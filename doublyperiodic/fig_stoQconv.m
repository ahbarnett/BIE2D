function fig_stoQconv
% Makes figs for doubly-periodic empty BVP convergence, Stokes.
% Barnett 5/28/16.   Adapted from fig_lapQconv.m. BIE2D 6/29/16

warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square (LU much worse)

mu = 0.7;   % overall viscosity const for Stokes PDE
rho = 1.5; z0 = rho*exp(0.7i); f0 = [0.3;-0.6]; % known Stokeslet loc & force

U.e1 = 1; U.e2 = 1i;    % unit cell lattice vectors
Rp = 1.4;               % proxy radius
M = 100; p.x = Rp * exp(1i*(0:M-1)'/M*2*pi); p = setupquad(p);    % proxy pts
proxyrep = @StoSLP;     % set proxy type via a kernel function

nt = 100; t.x = rand(nt,1)-0.5 + 1i*(rand(nt,1)-0.5);  % rand test pts in UC
%[xx yy] = meshgrid(linspace(-.5,.5,10)); t.x = xx(:)+1i*yy(:); % tp grid in UC
Atest = proxyrep(t,p,mu);   % evaluation matrix, only changes when M does

m = 22; [U L R B T] = doublywalls(U,m);
g = discrep(L,R,B,T,@(t) vTknown(t,z0,f0,mu));
w = [0*L.w 0*L.w L.w 0*L.w 0*B.w 0*B.w 0*B.w B.w]';   % 8m-cpt left null-vector?
nsings = 6;           % how many singular values to keep and show
Ms = 10:5:120;
verrs = nan*Ms'; sings = nan(numel(Ms),nsings); nrms = verrs; nwtqs = verrs;
for i=1:numel(Ms), M = Ms(i);        % ------ convergence in M (# proxy pts)
  p.x = Rp * exp(1i*(0:M-1)'/M*2*pi); p = setupquad(p);    % proxy pts
  Atest = proxyrep(t,p,mu);          % vel evaluation matrix
  Q = Qmat(p,L,R,B,T,proxyrep,mu);
  nwtqs(i) = norm(w'*Q);             % check presumptive left null vector
  S = svd(Q);
  sings(i,1) = max(S); sings(i,2:nsings) = S(end-nsings+2:end); % 1st & last few
  xi = linsolve(Q,g,lso); nrms(i) = norm(xi);
  verr = Atest*xi - vTknown(t,z0,f0,mu);        % vel err at all test pts
  verr(1:nt) = verr(1:nt)-verr(1);       % fix offset in x cpt from 1st test pt
  verr(nt+1:end) = verr(nt+1:end)-verr(1+nt);   % ..and y cpt
  verrs(i) = max(abs(verr));
  %if M==40, [U S V] = svd(Q); V(:,end), Atest*V(:,end), end  % sing vec is approx constant, and generates consts to high acc
end                                 % -------
figure; semilogy(Ms,[verrs],'+-'); hold on;
semilogy(Ms,sings(:,2:end),'-','color',.5*[1 1 1]);
semilogy(Ms,nrms,'go-');
% semilogy(Ms,nwtqs,'rs-'); % boring, always zero
rB = sqrt(0.5); semilogy(Ms,0.05*(rho/rB).^(-Ms/2),'r--');  % conv rate!
xlabel('M');
axis([min(Ms) max(Ms) 1e-17 max(nrms)]);
text(15,max(nrms)/10,'(c)','fontsize',12);
text(90,max(nrms)/10,sprintf('$m=%d$',m),'interpreter','latex','fontsize',12);
set(gcf,'paperposition',[0 0 3 3]);
set(gca,'ytickmode','manual', 'ytick',[1e-15 1e-10 1e-5 1 1e5]);
%print -depsc2 figs/stoQMconv.eps
%%%%%%%%%%%%%%%%%

function Q = Qmat(p,L,R,B,T,proxyrep,mu)  % matrix Q given proxy and colloc pts
[QL,~,QLn] = proxyrep(L,p,mu); [QR,~,QRn] = proxyrep(R,p,mu); % vel & trac, no p
[QB,~,QBn] = proxyrep(B,p,mu); [QT,~,QTn] = proxyrep(T,p,mu);
Q = [QR-QL; QRn-QLn; QT-QB; QTn-QBn];

function [ve Te] = vTknown(t,z0,f0,mu)  % known vel, trac on t, targ seg
if nargout==1, ve = StoSLP(t,struct('x',z0,'w',1),mu,f0);
else, [ve,~,Te] = StoSLP(t,struct('x',z0,'w',1),mu,f0);
end

function g = discrep(L,R,B,T,vTfun)   % discrepancy of given solution
% vTfun should have interface [velvec, tractionvec] = vTfun(targetsegment)
% Other inputs: LRBT are walls. Outputs: g is 8m col vec
[vR,tR] = vTfun(R); [vL,tL] = vTfun(L);
[vT,tT] = vTfun(T); [vB,tB] = vTfun(B);
g = [vR-vL; tR-tL; vT-vB; tT-tB];
