function fig_stoQconv
% Makes figs for doubly-periodic empty BVP convergence, Stokes.
% Barnett 5/28/16.   Adapted from fig_lapQconv.m. BIE2D 6/29/16

warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square (LU much worse)

mu = 0.7;   % overall viscosity const for Stokes PDE
rho = 1.5; z0 = rho*exp(0.7i); f0 = [0.3;-0.6]; % known Stokeslet loc & force

U.e1 = 1; U.e2 = 1i; U.nei = 0;  % unit cell
Rp = 1.4;
M = 100; p.x = Rp * exp(1i*(1:M)'/M*2*pi); p = quadr(p);    % proxy pts
proxyrep = @StoSLP;   % set proxy type via a kernel function

nt = 100; t.x = rand(nt,1)-0.5 + 1i*(rand(nt,1)-0.5);  % rand test pts in UC
%[xx yy] = meshgrid(linspace(-.5,.5,10)); t.x = xx(:)+1i*yy(:); % tp grid in UC
Atest = proxyrep(t,p,mu);   % evaluation matrix, only changes when M does

m = 20; [U L R B T] = doublywalls(U,m);
g = discrep(L,R,B,T,@(t) vTknown(t,z0,f0,mu));
w = [0*L.w 0*L.w L.w 0*L.w 0*B.w 0*B.w 0*B.w B.w]';   % 8m-cpt left null-vector?
nsings = 6;           % how many singular values to keep and show
Ms = 10:5:120;
verrs = nan*Ms'; sings = nan(numel(Ms),nsings); nrms = verrs; nwtqs = verrs;
for i=1:numel(Ms), M = Ms(i);        % ------ convergence in # proxy pts
  p.x = Rp * exp(1i*(1:M)'/M*2*pi); p = quadr(p); % proxy pts
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
text(90,max(nrms)/10,'$m=20$','interpreter','latex','fontsize',12);
set(gcf,'paperposition',[0 0 3 3]);
set(gca,'ytickmode','manual', 'ytick',[1e-15 1e-10 1e-5 1 1e5]);
%print -depsc2 ../paper/stoQMconv.eps
%%%%%%%%%%%%%%%%%

function Q = Qmat(p,L,R,B,T,proxyrep,mu)  % matrix Q given proxy and colloc pts
[QL QLn] = proxyrep(L,p,mu); [QR QRn] = proxyrep(R,p,mu);
[QB QBn] = proxyrep(B,p,mu); [QT QTn] = proxyrep(T,p,mu);
Q = [QR-QL; QRn-QLn; QT-QB; QTn-QBn];

function [ve Te] = vTknown(t,z0,f0,mu)  % known vel & traction on t targ seg
if nargout<2
  ve12 = SLPmatrix(t,struct('x',z0,'w',1),mu);
else
  [ve12 Te12] = SLPmatrix(t,struct('x',z0,'w',1),mu); Te = Te12*f0;
end
ve = ve12*f0;    % turn *-by-2 matrices into col vec using source force vec

function g = discrep(L,R,B,T,vTfun)   % discrepancy of given solution
% vTfun has interface [velvec, tractionvec] = vTfun(targetsegment)
% Other inputs: LRBT are walls. Outputs: g is 8m col vec
[vR tR] = vTfun(R); [vL tL] = vTfun(L);
[vT tT] = vTfun(T); [vB tB] = vTfun(B);
g = [vR-vL; tR-tL; vT-vB; tT-tB];

function s = quadr(s)  % interp curve geom, periodic trapezoid quadr
N = length(s.x); s.xp = perispecdiff(s.x); s.xpp = perispecdiff(s.xp);
s.sp = abs(s.xp); s.tang = s.xp./s.sp; s.nx = -1i*s.tang; 
s.cur = -real(conj(s.xpp).*s.nx)./s.sp.^2; s.w = 2*pi/N*s.sp; % speed weights
s.cw = 1i*s.nx.*s.w;  % complex weights (incl complex speed)

function g = perispecdiff(f)
% PERISPECDIFF - use FFT to take periodic spectral differentiation of vector
%
% g = PERISPECDIFF(f) returns g the derivative of the spectral interpolant
%  of f, which is assumed to be the values of a smooth 2pi-periodic function
%  at the N gridpoints 2.pi.j/N, for j=1,..,N (or any translation of such
%  points).
%
% Barnett 2/18/14
N = numel(f);
if mod(N,2)==0   % even
  g = ifft(fft(f(:)).*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].');
else
  g = ifft(fft(f(:)).*[0 1i*(1:(N-1)/2) 1i*((1-N)/2:-1)].');
end
g = reshape(g,size(f));

function h = showseg(s, U) % plot segments & possibly its 2d nei copies
if nargin<2, U.nei = 0; U.e1 = 0; U.e2 = 0; end % dummy uc
if iscell(s), for i=1:numel(s), showseg(s{i},U); end, return, end % multi segs
for i=-U.nei:U.nei, for j=-U.nei:U.nei
  plot(s.x+U.e1*i+U.e2*j, 'b.-'); axis equal xy;
  l=0.05; hold on; plot([s.x+U.e1*i+U.e2*j, s.x+l*s.nx+U.e1*i+U.e2*j].', 'k-');
end, end

function [A,T] = SLPmatrix(t,s,mu,a) % single-layer 2D Stokes kernel vel matrix
% Returns 2N-by-2N matrix from src force vector to 2 flow component
% t = target seg (x cols), s = src seg, a = optional translation of src seg
% No option for self-int, gives Inf on diag. Barnett 3/2/14
% 2nd output is traction matrix, needs t.nx normal (C-#); no self-int either.
if nargin<4, a = 0; end
N = numel(s.x); M = numel(t.x);
r = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
irr = 1./(conj(r).*r);    % 1/r^2
d1 = real(r); d2 = imag(r);
Ilogr = -log(abs(r));  % log(1/r) diag block
c = 1/(4*pi*mu);       % factor from Hsiao-Wendland book, Ladyzhenskaya
A12 = c*d1.*d2.*irr;   % off diag vel block
A = [c*(Ilogr + d1.^2.*irr), A12;                         % u_x
     A12,                c*(Ilogr + d2.^2.*irr)];         % u_y
A = A .* repmat([s.w(:)' s.w(:)'], [2*M 1]);              % quadr wei
if nargout>1           % traction (negative of DLP vel matrix w/ nx,ny swapped)
  rdotn = d1.*repmat(real(t.nx), [1 N]) + d2.*repmat(imag(t.nx), [1 N]);
  rdotnir4 = rdotn.*(irr.*irr); clear rdotn
  A12 = -(1/pi)*d1.*d2.*rdotnir4;
  T = [-(1/pi)*d1.^2.*rdotnir4,   A12;                   % own derivation
     A12,                      -(1/pi)*d2.^2.*rdotnir4];
  T = T .* repmat([s.w(:)' s.w(:)'], [2*M 1]);            % quadr wei
end

function [A,T] = DLPmatrix(t,s,mu,a) % double-layer 2D Stokes vel kernel matrix
% Returns 2N-by-2N matrix from src force vector to 2 flow components.
% If detects self-int, does correct diagonal limit, no jump condition (PV int).
% t = target seg (x cols), s = src seg, a = optional transl of src seg.
% 2nd output is optional traction matrix, without self-eval. Barnett 3/2/14
if nargin<4, a = 0; end
N = numel(s.x); M = numel(t.x);
r = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
irr = 1./(conj(r).*r);    % 1/R^2
d1 = real(r); d2 = imag(r);
rdotny = d1.*repmat(real(s.nx)', [M 1]) + d2.*repmat(imag(s.nx)', [M 1]);
rdotnir4 = rdotny.*(irr.*irr); if nargout<=1, clear rdotny; end
A12 = (1/pi)*d1.*d2.*rdotnir4;  % off diag vel block
A = [(1/pi)*d1.^2.*rdotnir4,   A12;                     % Ladyzhenskaya
     A12,                      (1/pi)*d2.^2.*rdotnir4];
if numel(s.x)==numel(t.x) && max(abs(s.x+a-t.x))<1e-14
  c = -s.cur/2/pi;           % diagonal limit of Laplace DLP
  tx = 1i*s.nx; t1=real(tx); t2=imag(tx);     % tangent vectors on src curve
  A(sub2ind(size(A),1:N,1:N)) = c.*t1.^2;       % overwrite diags of 4 blocks:
  A(sub2ind(size(A),1+N:2*N,1:N)) = c.*t1.*t2;
  A(sub2ind(size(A),1:N,1+N:2*N)) = c.*t1.*t2;
  A(sub2ind(size(A),1+N:2*N,1+N:2*N)) = c.*t2.^2;
end
A = A .* repmat([s.w(:)' s.w(:)'], [2*M 1]);            % quadr wei
if nargout>1           % traction, my formula
  nx1 = repmat(real(t.nx), [1 N]); nx2 = repmat(imag(t.nx), [1 N]);
  rdotnx = d1.*nx1 + d2.*nx2;
  ny1 = repmat(real(s.nx)', [M 1]); ny2 = repmat(imag(s.nx)', [M 1]);
  dx = rdotnx.*irr; dy = rdotny.*irr; dxdy = dx.*dy;
  R12 = d1.*d2.*irr; R = [d1.^2.*irr, R12; R12, d2.^2.*irr];
  nydotnx = nx1.*ny1 + nx2.*ny2;
  T = R.*kron(ones(2), nydotnx.*irr - 8*dxdy) + kron(eye(2), dxdy);
  T = T + [nx1.*ny1.*irr, nx1.*ny2.*irr; nx2.*ny1.*irr, nx2.*ny2.*irr] + ...
      kron(ones(2),dx.*irr) .* [ny1.*d1, ny1.*d2; ny2.*d1, ny2.*d2] + ...
      kron(ones(2),dy.*irr) .* [d1.*nx1, d1.*nx2; d2.*nx1, d2.*nx2];
  T = (mu/pi) * T;                                        % prefac
  T = T .* repmat([s.w(:)' s.w(:)'], [2*M 1]);            % quadr wei
end

function i = diagind(A) % return indices of diagonal of square matrix
N = size(A,1); i = sub2ind(size(A), 1:N, 1:N);

% GAUSS  nodes x (Legendre points) and weights w
%        for Gauss quadrature on [-1,1], for N small (<100). Trefethen book.
function [x,w] = gauss(N)
beta = .5./sqrt(1-(2*(1:N-1)).^(-2));
T = diag(beta,1) + diag(beta,-1);
[V,D] = eig(T);
x = diag(D); [x,i] = sort(x);
w = 2*V(1,i).^2;
