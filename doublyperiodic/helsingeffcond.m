  function effcond = helsingeffcond(c,sigma2,npan0,verb)
% HELSINGEFFCOND  Solve square disc lattice effective conductivity using RCIP.
%
% effcond = helsingeffcond(c,sigma2,npan0)
% effcond = helsingeffcond(c,sigma2,npan0,verb)
%
%  Uses demo14b code from Helsing's RCIP tutorial to return effcond for
%  discs of radius r in a 1-periodic array in the plane. Only half the number
%  of refinements as Helsing are used, since theory gives that the densities
%  are smooth on the scale 1/c, which is much larger than 1/c^2.
%
% Inputs:
%  c - radius parameter (inverse dist to singularity): r = sqrt(1-1/c^2)/2.
%  sigma2 - conductivity sigma inside disc (exterior sigma is unity).
%  npan0 (optional) - number of panels per 1/4 circle; default 4
%                     (Helsing checks that already converged by this value).
%  verb (optional) - 0:silent, 1:text output.
% Output:
%  effcond - effective conductivity. Uses same applied gradient as demo14b.m
%
% Barnett 9/27/17
% Simply wrapping & tweaking Helsing's cool demo14b.m code of 5/19/17

% **************************************************************
% *** Last changed 2017-05-19 by Johan Helsing                 *
% *** Effective conductivity of doubly periodic array of disks *
% **************************************************************

if nargin<3 || isempty(npan0), npan0 = 4; end
if nargin<4, verb = 0; end
  Ng=22;           % number of nodes on canonical interval (ie, pan order, ahb)
  [W,T]=gaussP(Ng);
  [IP,IPW]=IPinit(T,W);
  Pbc=Pbcinit(IP);
  PWbc=Pbcinit(IPW);
%  evec=1;                      % applied electric field
%  evec=1i                      % applied electric field
  evec=(1+1i)/sqrt(2);         % applied electric field  % why? isotropic
  sigma1=1;                    % conductivity of background medium  
  lambda=(sigma2-sigma1)./(sigma2+sigma1);   % good
  d=1./(c.^2+c.*sqrt(c.^2-1));      % closest distance
  rr=sqrt(c.^2-1)./(2*c);              % radius (for unit periodicity)
  if verb, fprintf('radius = %.15g\n',rr), end
  %effcond_asymp=McPhedran(c,lambda,sigma1);
  nplot=1;   % hack
  relerr=zeros(nplot,3);
  npvec=zeros(nplot,1);
  for k=1:nplot               
    %npan0=npan0vec(k);         % number of panels on 1/4 circle
    np0=Ng*npan0;
    parlen=2/(4*npan0)*pi;     % panel length measured in parameter
    npan=4*npan0;              % number of panels on circle
    if verb, fprintf('N coarse nodes on whole circle = %g\n',npan*Ng), end % ahb
    sinter=pi*linspace(-1,1,npan+1)';
    sinterdiff=parlen*ones(npan,1);
    izone=zeros(4*Ng,4);
    izone(:,1)=[4*np0-2*Ng+1:4*np0 1:2*Ng]';
    izone(:,2)=(  np0-2*Ng+1:  np0+2*Ng)';
    izone(:,3)=(2*np0-2*Ng+1:2*np0+2*Ng)';
    izone(:,4)=(3*np0-2*Ng+1:3*np0+2*Ng)';
    for kk=1 %1:3       % just do one expt - hack
      panlen=rr(kk)*parlen;        % actual panel length on coarse mesh
      nsub=ceil(.5*log(panlen/d(kk))/log(2));  % Johan chose how many subdivisions based on d scale, but ahb halved it: density changes only on sqrt(d)~1/c scale.
      if nsub<1, nsub=1; end        % handle d being too large
      if verb, fprintf('nsub = %d\n',nsub), end
      [z,zp,zpp,w,wzp,np,zinter]=zinit(sinter,sinterdiff,T,W,rr(kk),npan);
      K=MRinitper(z,zp,zpp,w,wzp,zinter,izone,npan0,Ng,np);
      R=eye(np);
      Rblock=Rcomp(T,W,Pbc,PWbc,parlen,lambda(kk),rr(kk),nsub,d(kk),Ng);
      myind1=1:4*Ng;
      myind2=4*Ng+1:8*Ng;
      for ic=1:4
	R(izone(:,ic),izone(:,ic))=Rblock(myind1,myind1);
      end
      R(izone(:,1),izone(:,3))=Rblock(myind2,myind1);
      R(izone(:,3),izone(:,1))=Rblock(myind1,myind2);
      R(izone(:,2),izone(:,4))=Rblock(myind2,myind1);      
      R(izone(:,4),izone(:,2))=Rblock(myind1,myind2);
      R=sparse(R);
      g=2*lambda(kk)*imag(conj(evec)*z);     
      [rhotilde,it]=myGMRESR(K,R,lambda(kk),g,np,60,eps,verb);
      if verb, disp(['GMRES_iter=',num2str(it)]), end
      rhohat=R*rhotilde;
      effcond(kk)=sigma1*(1-rhohat.'*real(conj(evec)*wzp));
      %relerr(k,kk)=abs(effcond(kk)-effref(kk))/abs(effref(kk)); % is no ref
    end
    npvec(k)=np;
    %myplot(relerr,npvec,k)  
    %effcond is returned
  end

  function sigeff=McPhedran(c,lambda,sigma1)
% **************************************************************
% *** Asymptotic formula from:                                 *
% *** R.C. McPhedran, L. Poladian, G.W. Milton, `Asymptotic    *
% *** studies of closely spaced, highly conducting cylinders', *
% *** Proc. R. Soc. Lond. A{\bf 415}, 185--196 (1988).         *
% **************************************************************
  eulmasc=0.577215664901532860606512090082402431042;
  s=log(lambda)./log((c-1)./(c+1));
  sigeff=sigma1*(1+pi*(c-1)./(1+2*s.*(log(c)-eulmasc-psi(1+s))));
  
  function [z,zp,zpp,w,wzp,zinter]=zlocinit6(T,W,rr,sinterdiff)
  w=[W/4;W/4;W/2]*sinterdiff;
  w=[flipud(w);w];
  s=[(T+1)/4;(T+3)/4;(T+3)/2]*sinterdiff;
  s=[-flipud(s);s];
  z=rr*(-2*sin(s/2).^2+1i*sin(s));
  zp=1i*rr*exp(1i*s);
  zpp=-rr*exp(1i*s);
  wzp=w.*zp;
  sinter=[-2;-1;-0.5;0;0.5;1;2]*sinterdiff;
  zinter=rr*(-2*sin(sinter/2).^2+1i*sin(sinter));

  function [z,zp,zpp,w,wzp,np,zinter]=zinit(sinter,sinterdiff,T,W,rr,npan)
  Ng=length(T);
  np=Ng*npan;
  s=zeros(np,1);
  w=zeros(np,1);
  for k=1:npan
    myind=(k-1)*Ng+1:k*Ng;
    sdif=sinterdiff(k)/2;
    s(myind)=(sinter(k)+sinter(k+1))/2+sdif*T;
    w(myind)=W*sdif;
  end
  z=zfunc(s,rr);
  zp=zpfunc(s,rr);
  zpp=zppfunc(s,rr);
  wzp=w.*zp;
  zinter=zfunc(sinter,rr);  

  function zout=zfunc(s,rr)
  zout=rr*exp(1i*s);

  function zpout=zpfunc(s,rr)
  zpout=1i*rr*exp(1i*s);

  function zppout=zppfunc(s,rr)
  zppout=-rr*exp(1i*s);
  
  function R12=Rcomp(T,W,Pbc,PWbc,sinterdiff0,lambda,rr,nsub,d,Ng)
  starL=Ng+1:5*Ng;
  circL=[1:Ng 5*Ng+1:6*Ng];
  starS=Ng+1:3*Ng;
  circS=[1:Ng 3*Ng+1:4*Ng];
  starL12=[starL 6*Ng+starL];
  circL12=[circL 6*Ng+circL];
  starS12=[starS 4*Ng+starS];
  circS12=[circS 4*Ng+circS];
  Pbc12=blkdiag(Pbc,Pbc);
  PWbc12=blkdiag(PWbc,PWbc);
  for level=1:nsub
    sinterdiff=sinterdiff0/2^(nsub-level);
    [z,zp,zpp,w,wzp,zinter]=zlocinit6(T,W,rr,sinterdiff);
    K=MRinit(z,zp,zpp,w,wzp,6*Ng);
    z1= z-d/2;
    z2=-z+d/2;
    wzp1= wzp;
    wzp2=-wzp;
    zinter1= zinter-d/2;
    zinter2=-zinter+d/2;
    K12=zeros(12*Ng);
    K12(     1: 6*Ng,     1: 6*Ng)=K;
    K12(6*Ng+1:12*Ng,6*Ng+1:12*Ng)=K;    
    if level==1
      K12(     1: 6*Ng,6*Ng+1:12*Ng)=MRinit12(z2,z1,zinter2,wzp2,6,Ng);
      K12(6*Ng+1:12*Ng,     1: 6*Ng)=MRinit12(z1,z2,zinter1,wzp1,6,Ng);
      MAT12=eye(12*Ng)+lambda*K12;   
      R12=inv(MAT12(starL12,starL12));
    else
      K12(     1: 6*Ng,6*Ng+1:12*Ng)=MRinit12s(z2,z1,zinter2,wzp2,Ng);
      K12(6*Ng+1:12*Ng,     1: 6*Ng)=MRinit12s(z1,z2,zinter1,wzp1,Ng);
      MAT12=eye(12*Ng)+lambda*K12;   
    end
    R12=SchurBana(Pbc12,PWbc12,MAT12,R12,starL12,circL12,starS12,circS12);
  end

  function Pbc=Pbcinit(IP)
  Pbc=blkdiag(IP,rot90(IP,2));
  
  function A=SchurBana(P,PW,K,A,starL,circL,starS,circS)
  VA=K(circL,starL)*A;
  PTA=PW'*A;
  PTAU=PTA*K(starL,circL);
  DVAUI=inv(K(circL,circL)-VA*K(starL,circL));
  DVAUIVAP=DVAUI*(VA*P);
  A(starS,starS)=PTA*P+PTAU*DVAUIVAP;
  A(circS,circS)=DVAUI;
  A(circS,starS)=-DVAUIVAP;
  A(starS,circS)=-PTAU*DVAUI;

  function M1=MRinit(z,zp,zpp,w,wzp,N)
% *** double layer potential (Neumann-Poincaré operator) ***
  M1=zeros(N);
  for m=1:N
     M1(:,m)=imag(wzp(m)./(z(m)-z));
  end
  M1(1:N+1:N^2)=w.*imag(zpp./zp)/2;
  M1=M1/pi;
  
  function M1=MRinitper(z,zp,zpp,w,wzp,zinter,izone,npan0,Ng,N)
  ztrans=[1+1i;1;1-1i;1i;-1i;-1+1i;-1;-1-1i;0];
  M1=MRinit(z,zp,zpp,w,wzp,N);
  for k=1:4
    M1(izone(:,k),izone(:,k))=0;
  end
  for k=1:8
    if k==1 || k==3 || k==6 || k==8
      TMP=zeros(N);
      for n=1:N
	TMP(:,n)=TMP(:,n)+imag(wzp(n)./(z(n)+ztrans(k)-z))/pi;
      end
    elseif k==2
      npans=[4*npan0 1 2*npan0 2*npan0+1];
      TMP=MRinit12p(z+ztrans(k),z,zinter+ztrans(k),wzp,4*npan0,Ng,npans);
      TMP(izone(:,3),izone(:,1))=0;
    elseif k==4
      npans=[npan0 npan0+1 3*npan0 3*npan0+1];
      TMP=MRinit12p(z+ztrans(k),z,zinter+ztrans(k),wzp,4*npan0,Ng,npans);
      TMP(izone(:,4),izone(:,2))=0;      
    elseif k==5
      npans=[npan0 npan0+1 3*npan0 3*npan0+1];
      TMP=MRinit12p(z+ztrans(k),z,zinter+ztrans(k),wzp,4*npan0,Ng,npans);
      TMP(izone(:,2),izone(:,4))=0;      
    elseif k==7
      npans=[4*npan0 1 2*npan0 2*npan0+1];
      TMP=MRinit12p(z+ztrans(k),z,zinter+ztrans(k),wzp,4*npan0,Ng,npans);
      TMP(izone(:,1),izone(:,3))=0;      
    end
    M1=M1+TMP;
  end
  SS=SSinit;
  for n=1:N
    zd=z(n)-z;
    zd4=zd.^4;
    zsum=SS(104);
    for k=25:-1:1
      zsum=SS(4*k)+zsum.*zd4;
    end
    zsum=zsum.*zd.^3-conj(zd)*SS(2);
    M1(:,n)=M1(:,n)-imag(wzp(n)*zsum)/pi;  
  end
  
  function M1=MRinit12(zs,zt,zinters,wzps,npan,Ng,npan0)
  icv=ones(npan,1);
  if nargin==7
    icv([npan0:npan0+1 3*npan0:3*npan0+1])=0;  
  end
  M1=zeros(Ng*npan);
  for n=1:Ng*npan
     M1(:,n)=imag(wzps(n)./(zs(n)-zt));
  end
  for n1=1:npan
    n=(n1-1)*Ng+(1:Ng);
    b=zinters(n1);
    c=zinters(n1+1);
    cmb=c-b;
    cpb=c+b;
    rjtr=(2*zs(n)-cpb)/cmb;
    A=ones(Ng);
    for kx=2:Ng
      A(:,kx)=rjtr.*A(:,kx-1);
    end
    for m1=1:npan
      for m2=1:Ng
	m=(m1-1)*Ng+m2;
	rtr=(2*zt(m)-cpb)/cmb;
	if abs(rtr)<2.0 && icv(n1)==1 && icv(m1)==1
	  M1(m,n)=imag(wCinit(A,rtr,Ng));
	end
      end
    end
  end
  M1=M1/pi;
  
  function M1=MRinit12p(zs,zt,zinters,wzps,npan,Ng,npans)
  icv=ones(npan,1);
  icv(npans)=0;  
  M1=zeros(Ng*npan);
  for n=1:Ng*npan
     M1(:,n)=imag(wzps(n)./(zs(n)-zt));
  end
  for n1=1:npan
    n=(n1-1)*Ng+(1:Ng);
    b=zinters(n1);
    c=zinters(n1+1);
    cmb=c-b;
    cpb=c+b;
    rjtr=(2*zs(n)-cpb)/cmb;
    A=ones(Ng);
    for kx=2:Ng
      A(:,kx)=rjtr.*A(:,kx-1);
    end
    for m1=1:npan
      for m2=1:Ng
	m=(m1-1)*Ng+m2;
	rtr=(2*zt(m)-cpb)/cmb;
	if abs(rtr)<2.0 && icv(n1)==1 && icv(m1)==1
	  M1(m,n)=imag(wCinit(A,rtr,Ng));
	end
      end
    end
  end
  M1=M1/pi;
  
  function M1=MRinit12s(zs,zt,zinters,wzps,Ng)
  IC=ones(6);
  IC(2:5,2:5)=0;
  M1=zeros(Ng*6);
  for n=1:Ng*6
     M1(:,n)=imag(wzps(n)./(zs(n)-zt));
  end
  for n1=[1 2 5 6]
    n=(n1-1)*Ng+(1:Ng);
    b=zinters(n1);
    c=zinters(n1+1);
    cmb=c-b;
    cpb=c+b;
    rjtr=(2*zs(n)-cpb)/cmb;
    A=ones(Ng);
    for kx=2:Ng
      A(:,kx)=rjtr.*A(:,kx-1);
    end
    for m1=[1 2 5 6]
      if IC(n1,m1)==1
	for m2=1:Ng
	  m=(m1-1)*Ng+m2;
	  rtr=(2*zt(m)-cpb)/cmb;
	  if abs(rtr)<2.0
	    M1(m,n)=imag(wCinit(A,rtr,Ng));
	  end
	end
      end
    end
  end
  M1=M1/pi;
  
  function M2=wCinit(A,rtr,Ng)
  p=zeros(Ng+1,1);
  c=(1-(-1).^(1:Ng))./(1:Ng);
  upp=log( 1-rtr);
  loo=log(-1-rtr);
  if imag(rtr)>0 && abs(real(rtr))<1
    loo=loo+2i*pi;
  end
  p(1)=upp-loo;
  for k=1:Ng
    p(k+1)=p(k)*rtr+c(k);	% forward recursion
  end
  M2=A.'\p(1:Ng);

  function [x,it]=myGMRESR(A,R,lmb,b,n,m,tol,verb)
% *** GMRES with low-threshold stagnation control ***
  V=zeros(n,m+1);
  H=zeros(m);
  cs=zeros(m,1);
  sn=zeros(m,1);
  bnrm2=norm(b);
  V(:,1)=b/bnrm2;
  s=bnrm2*eye(m+1,1);
  for it=1:m                                  
    it1=it+1;                                   
    w=lmb*(A*(R*V(:,it)));
    for k=1:it
      H(k,it)=V(:,k)'*w;
      w=w-H(k,it)*V(:,k);
    end
    H(it,it)=H(it,it)+1;
    wnrm2=norm(w);
    V(:,it1)=w/wnrm2;
    for k=1:it-1                                
      temp     = cs(k)*H(k,it)+sn(k)*H(k+1,it);
      H(k+1,it)=-sn(k)*H(k,it)+cs(k)*H(k+1,it);
      H(k,it)  = temp;
    end
    [cs(it),sn(it)]=rotmat(H(it,it),wnrm2);     
    H(it,it)= cs(it)*H(it,it)+sn(it)*wnrm2;
    s(it1) =-sn(it)*s(it);                      
    s(it)  = cs(it)*s(it);                         
    myerr=abs(s(it1))/bnrm2;
    if (myerr<=tol)||(it==m)                     
      if verb, disp(['predicted residual = ' num2str(myerr)]), end
      y=triu(H(1:it,1:it))\s(1:it);             
      x=fliplr(V(:,1:it))*flipud(y);
      trueres=norm(x+lmb*(A*(R*x))-b)/bnrm2;
      if verb, disp(['true residual      = ',num2str(trueres)]), end
      break
    end
  end
  
  function [c,s]=rotmat(a,b)
% *** Compute the Givens rotation matrix parameters for a and b ***
  if  b==0
    c=1;
    s=0;
  elseif abs(b)>abs(a)
    temp=a/b;
    s=1/sqrt(1+temp^2);
    c=temp*s;
  else
    temp=b/a;
    c=1/sqrt(1+temp^2);
    s=temp*c;
  end

  function [IP,IPW]=IPinit(T,W)
  Ng=length(T);
  A=ones(Ng);
  AA=ones(2*Ng,Ng);
  T2=[T-1;T+1]/2;
  W2=[W;W]/2;
  for k=2:Ng
    A(:,k)=A(:,k-1).*T;
    AA(:,k)=AA(:,k-1).*T2;   
  end
  IP=AA/A;
  IPW=IP.*(W2*(1./W)');

  function myplot(relerr,npvec,np)
  relerr(relerr<eps)=eps;
  figure(1)
  loglog(npvec(1:np),relerr(1:np,1),'bo',npvec(1:np),relerr(1:np,2),'rv', ...
	 npvec(1:np),relerr(1:np,3),'g*')
  hl=legend('$\sigma_2=10^8$; $c=3\cdot 10^3$', ...
	    '$\sigma_2=10^4$; $c=10^4$', ...
	    '$\sigma_2=10^3$; $c=10^3$');
  set(hl,'Interpreter','LaTeX','FontSize',12)    
  grid on
  title('Convergence of $\sigma_{\rm eff}$ with mesh refinement', ...
	'Interpreter','LaTeX')
  ylabel('estimated relative error in $\sigma_{\rm eff}$', ...
	 'Interpreter','LaTeX')
  xlabel('number of discretization points on coarse grid', ...
	 'Interpreter','LaTeX')
  axis([3e2 1e4 1e-16 1e-0])
  axis square
  drawnow
    
  function [W,T]=gaussP(N)
  k=1:N-1;
  B=k./sqrt(4*k.*k-1);
  [V,D]=eig(diag(B,1)+diag(B,-1));
  DD=diag(D);
  [~,myind0]=sort(DD);
  T=DD(myind0);
  W=2*(V(1,myind0).*V(1,myind0))';
  
  function SS=SSinit
  SS=zeros(104,1);
  SS(2)  =3.1415926535897932384626433832795028;
  SS(4)  =0.15121200215389753821768994224868884;
  SS(8)  =5.7730353651895184471546807401316399e-03;
  SS(12) =1.3490128279703747512623657685445488e-03;
  SS(16) =7.0033025024855874182489529648906547e-05;
  SS(20) =3.0031762895595772422633093074649335e-06;
  SS(24) =2.4280383862810120496714860850943050e-07;
  SS(28) =1.6099409965709787749402945332222719e-08;
  SS(32) =8.9758066666893286269089449457352552e-10;
  SS(36) =5.7045151271059136564668362798495565e-11;
  SS(40) =3.7180310586178652778795231479449756e-12;
  SS(44) =2.2744025145537515649529231310184412e-13;
  SS(48) =1.4081285775545639932256139922281536e-14;
  SS(52) =8.9097416933543890760423470380997271e-16;
  SS(56) =5.5655838097435282723971282718107241e-17;
  SS(60) =3.4617327481687556094172723767053656e-18;
  SS(64) =2.1678173400149605207446309923699580e-19;
  SS(68) =1.3566184768600140463260467809434722e-20;
  SS(72) =8.4682093760432185521906660944906355e-22;
  SS(76) =5.2922456038352385260388276910920827e-23;
  SS(80) =3.3094447765677773871482964155303854e-24;
  SS(84) =2.0680633809167300063750942769111378e-25;
  SS(88) =1.2923290806020022572336690888375830e-26;
  SS(92) =8.0780717135949195591966820055206958e-28;
  SS(96) =5.0489043219484377842248498898718264e-29;
  SS(100)=3.1553782794313412460626512625194132e-30;
  SS(104)=1.9721357749725043881844218286505179e-31;
