function fig_lapQconv
% make figs for doubly periodic empty BVP convergence, Laplace case.
% All matrix-filling, native quadr.
% Barnett 5/9/16, adapted from perineu2dnei1.m. bie2d 6/29/16

warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square (LU much worse)

rho=1.5; z0 = rho * exp(0.7i);       % singularity in v exact soln
ve = @(z) log(abs(z-z0));   % exact soln and its partials...
vex = @(z) real(z-z0)./abs(z-z0).^2; vey = @(z) imag(z-z0)./abs(z-z0).^2;
disp('test partials of known exact v:'); testpartials(ve,vex,vey)

U.e1 = 1; U.e2 = 1i;      % unit cell lattice vectors
Rp = 1.4;                 % proxy radius
M = 100; p.x = Rp * exp(1i*(0:M-1)'/M*2*pi); p = setupquad(p);    % proxy pts
proxyrep = @LapSLP;       % sets proxy pt type via a kernel function call

nt = 100; t.x = rand(nt,1)-0.5 + 1i*(rand(nt,1)-0.5);  % rand test pts in UC
%[xx yy] = meshgrid(linspace(-.5,.5,10)); t.x = xx(:)+1i*yy(:); % tp grid in UC
Atest = proxyrep(t,p);   % evaluation matrix, only changes when M does

% m-conv plot (# wall pts)
ms = 2:2:24; verrs = nan*ms;
for i=1:numel(ms)
  [U L R B T] = doublywalls(U,ms(i));
  Q = Qmat(p,L,R,B,T,proxyrep);
  g = discrep(L,R,B,T,ve,vex,vey);
  xi = linsolve(Q,g,lso); %norm(xi)
  verr = Atest*xi - ve(t.x);  % err at all test pts
  verr = verr-verr(1);   % fix the offset using the first test pt
  verrs(i) = max(abs(verr));
end
%figure; showsegment({p L R T B}); hold on; plot(z0,'*');
figure; semilogy(ms,verrs,'+-'); xlabel('m'); ylabel('max abs error in v');
axis([min(ms) max(ms) 1e-16 1]);
text(4,1e-1,'(a)','fontsize',12);
text(15,1e-1,sprintf('$M=%d$',M),'interpreter','latex','fontsize',12);
set(gcf,'paperposition',[0 0 3 3]);
%print -depsc2 figs/lapQmconv.eps

m = 22; [U L R B T] = doublywalls(U,m); g = discrep(L,R,B,T,ve,vex,vey); % fix
w = [0*L.w L.w 0*B.w B.w]';   % supposed left null-vector (note v will not be)
nsings = 6;           % how many singular values to keep and show
Ms = 10:5:120;
verrs = nan*Ms'; sings = nan(numel(Ms),nsings); nrms = verrs; nwtqs = verrs;
for i=1:numel(Ms), M = Ms(i);        % # proxy pts conv
  p.x = Rp * exp(1i*(1:M)'/M*2*pi); p = setupquad(p); % proxy pts
  Atest = proxyrep(t,p);   % evaluation matrix
  Q = Qmat(p,L,R,B,T,proxyrep);
  nwtqs(i) = norm(w'*Q);        % check left null vector
  S = svd(Q);
  %  nkeep = min(numel(S),nsings); sings(i,1:nkeep) = S(end-nkeep+1:end);
  sings(i,1) = max(S); sings(i,2:nsings) = S(end-nsings+2:end); % 1st & last few
  xi = linsolve(Q,g,lso); nrms(i) = norm(xi);
  verr = Atest*xi - ve(t.x);  % err at all test pts
  verr = verr-verr(1);   % fix the offset using the first test pt
  verrs(i) = max(abs(verr));
  %if M==40, [U S V] = svd(Q); V(:,end), Atest*V(:,end), end  % sing vec is approx constant, and generates consts to high acc
end
figure; semilogy(Ms,[verrs],'+-'); hold on;
semilogy(Ms,sings(:,2:end),'-','color',.5*[1 1 1]);
semilogy(Ms,nrms,'go-');
% semilogy(Ms,nwtqs,'rs-'); % boring, always zero
rB = sqrt(0.5); semilogy(Ms,0.05*(rho/rB).^(-Ms/2),'r--');  % conv rate!
xlabel('M');
%ylabel('max abs error in v');
axis([min(Ms) max(Ms) 1e-17 max(nrms)]);
text(15,max(nrms)/10,'(b)','fontsize',12);
text(90,max(nrms)/10,sprintf('$m=%d$',m),'interpreter','latex','fontsize',12);
set(gcf,'paperposition',[0 0 3 3]);
set(gca,'ytickmode','manual', 'ytick',[1e-15 1e-10 1e-5 1 1e5]);
%print -depsc2 figs/lapQMconv.eps

  
%%%%%%%%%%%%%%%%%

function testpartials(v,vx,vy)  % finite differencing test of analytic partials
eps = 1e-5;
z = 1-0.5i;
vx(z) - (v(z+eps)-v(z-eps))/(2*eps)        % should be around 1e-10
vy(z) - (v(z+1i*eps)-v(z-1i*eps))/(2*eps)  % "

function Q = Qmat(p,L,R,B,T,proxyrep) % matrix Q given proxy and colloc pts
[QL QLn] = proxyrep(L,p); [QR QRn] = proxyrep(R,p);
[QB QBn] = proxyrep(B,p); [QT QTn] = proxyrep(T,p);
Q = [QR-QL; QRn-QLn; QT-QB; QTn-QBn];

function g = discrep(L,R,B,T,v,vx,vy)   % discrepancy of known solution
% v, vx, vy are funcs from C to R. Output is 4m col vec
g1 = v(R.x)-v(L.x);
g2 = vx(R.x)-vx(L.x); % assumes vert
g3 = v(T.x)-v(B.x);
g4 = vy(T.x)-vy(B.x);  % assumes horiz
g = [g1;g2;g3;g4];

