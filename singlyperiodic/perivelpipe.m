function perivelpipe
% Longitudinal periodize 2D velocity-BC Stokes "pipe" geom w/ press drop (pgro)
% using circle of SLP proxy sources. Barnett, for Veerapaneni & Gillman 3/11/14

clear; v=2; expt='t'; % verbosity=0,1,2. expt='t' test, 'd' driven no-slip demo
makefigs = 0;   % whether to write to EPS files for writeup
uc.d = 2*pi; % unitcell, d=period in space
N = 100; %80; % pts per top and bottom wall (enough for 1e-15 in this case)
U.Z = @(t) t + 1i*(1.5+sin(t)); U.Zp = @(t) 1 + 1i*cos(t);
U.Zpp = @(t) -1i*sin(t); U = quadr(U,N); % t=0 is left end, t=2pi right end
D.Z = @(t) t + 1i*(-1.5+cos(2*t)); D.Zp = @(t) 1 - 2i*sin(2*t);
D.Zpp = @(t) - 4i*cos(2*t); D = quadr(D,N); % same direction (left to right)
U.nx = -U.nx; U.cur = -U.cur; % correct for sense of U, opp from periodicdirpipe
% normals on pipe walls point outwards
inside = @(z) imag(z-D.Z(real(z)))>0 & imag(z-U.Z(real(z)))<0; % ok U,D graphs

mu = 0.7;                                           % viscosity
if expt=='t' % Exact soln: either periodic or plus fixed pressure drop / period:
  %ue = @(x) [1+0*x; -2+0*x]; pe = @(x) 0*x; % the exact soln: uniform rigid flow, constant pressure everywhere (no drop)
  h=.2; ue = @(x) h*[imag(x).^2;0*x]; pe = @(x) h*2*mu*real(x); % horiz Poisseuil flow (pres drop)
  rhs = [ue(U.x); ue(D.x)]; % Driving: U,D bdry vels, each stacked as [u_1;u_2]
  pgro = pe(uc.d)-pe(0); % known pressure growth across one period (a number)
elseif expt=='d'
  rhs = zeros(4*N,1);       % no-slip BCs, ie homog vel data on U,D
  pgro = -1;                 % given pressure driving, for flow +x (a number)
end

uc.nei = 1; % how many nei copies either side (use eg 1e3 to test A w/o AP)
n = 30; % pts per side
[x w] = gauss(n); x = (1+x)/2; w = w'/2; % quadr on [0,1]
H = U.Z(0)-D.Z(0); L.x = D.Z(0) + H*x; L.nx = 0*L.x+1; L.w = H*w; % left side
R = L; R.x = L.x+uc.d; % right side
M = 50; % # proxy pts in periodizing basis (2 force comps per pt, so 2M dofs)
b.x = pi + 1.1*2*pi*exp(2i*pi*(1:M)'/M);        % the proxy pts
%b.nx = exp(2i*pi*(1:M)'/M);                % only needed if DLP proxy basis
b.w = 1+0*b.x;  % unit quadr weights (dummy)

nx = 60; gx = 2*pi*((1:nx)-0.5)/nx; ny = nx; gy = gx - pi; % plotting grid
[xx yy] = meshgrid(gx,gy); t.x = xx(:)+1i*yy(:); Mt = numel(t.x);
di = reshape(inside(t.x),size(xx));  % boolean if inside domain
if expt=='t', ueg = ue(t.x);             % evaluate exact soln on grid
  ue1 = reshape(ueg(1:Mt),size(xx)); ue2 = reshape(ueg(Mt+(1:Mt)),size(xx));
  peg = reshape(pe(t.x),size(xx));
  if v, figure; imagesc(gx,gy, peg); colormap(jet(256)); colorbar;
    hold on; quiver(gx,gy, ue1,ue2, 10); end
else, figure; end
if v, showseg(U,uc); hold on; showseg(D,uc); showseg(L); showseg(R);
  vline([0 uc.d]); axis xy equal tight;
  text(-.5,0,'L');text(2*pi+.5,0,'R');text(4,-1,'D');text(pi,0.5,'U');
  plot(b.x, 'r.'); title('geom');
  if expt=='t',title('geom and (u,p) known soln'); end
  if makefigs, axis([-4 10 -7 7]); print -depsc2 geom.eps, end
end

% fill system matrices A B C Q (subblock ordering always U then D), same
A = -eye(4*N)/2; % A's jump relation part. A maps density to vel vec, on U,D
for i=-uc.nei:uc.nei, a = i*uc.d;
  A = A + [DLPmatrix(U,U,mu,a) DLPmatrix(U,D,mu,a); DLPmatrix(D,U,mu,a) DLPmatrix(D,D,mu,a)]; end %figure; imagesc(A); colorbar; title('A');
% single layer (monopoles) on proxy points:
B = [SLPmatrix(U,b,mu); SLPmatrix(D,b,mu)]; % maps 2M peri dofs to vels on U,D
a=uc.nei*uc.d; [RU RUn] = DLPmatrix(R,U,mu,-a); [LU LUn] = DLPmatrix(L,U,mu,a);
[RD RDn] = DLPmatrix(R,D,mu,-a); [LD LDn] = DLPmatrix(L,D,mu,a);
C = [RU-LU, RD-LD; RUn-LUn, RDn-LDn]; % maps cancelled densities to discrepancy
[Rb Rbn] = SLPmatrix(R,b,mu); [Lb Lbn] = SLPmatrix(L,b,mu);
Q = [Rb-Lb; Rbn-Lbn]; % maps periodizing dofs to discrepancy

QdagC = Q\C;      % backwards stable least-sq solve
AP = A - B*QdagC; % system mat maps periodized density on U,D to vels on U,D
if v>1, figure; imagesc(AP); colorbar; title('A_P'); end

Tgro = -pgro * [real(R.nx);imag(R.nx)]; % traction driving growth (vector func)
Qdagg = Q\[zeros(2*n,1); Tgro]; % no vel growth; here []=g is rhs for peri rows
rhsgro = -B*Qdagg;       % change in rhs due to peri (along-pipe) driving
rhs = rhs + rhsgro;

% dense solve, note nullity(AP)=1 corresp to unknown pres const, but consistent:
tau = AP\rhs;          % tau is DLP "vector density" (ie bdry force dipole)
c = -QdagC*tau + Qdagg;        % get periodizing dofs (incl peri driving part)

if v, figure; plot([tau;c],'+-'); legend('[\tau;c]'); title('soln density and periodizing dofs'); end

if expt=='t', z = []; z.x = 2+1i;    % pointwise test u soln (target object z)
  u = SLPmatrix(z,b,mu) * c; % eval @ test pt: first do MFS (proxy) contrib
  for i=-uc.nei:uc.nei, a = i*uc.d;  % add in 3 copies of LPs on U,D...
    u = u + [DLPmatrix(z,U,mu,a) DLPmatrix(z,D,mu,a)] * tau;
  end
  fprintf('u error at pt = %.3g, ||tau||=%.3g\n', norm(u-ue(z.x)),norm(tau))
end

ug = SLPmatrix(t,b,mu) * c;    % eval on grid: MFS (proxy) contrib
pg = SLPpresmatrix(t,b,mu) * c;  % (vel and pres parts separate matrices)
for i=-uc.nei:uc.nei, a = i*uc.d;  % add in 3 copies of LPs on U,D...
  ug = ug + [DLPmatrix(t,U,mu,a) DLPmatrix(t,D,mu,a)] * tau; % uses RAM
  pg = pg + [DLPpresmatrix(t,U,mu,a) DLPpresmatrix(t,D,mu,a)] * tau;
end
u1 = reshape(ug(1:Mt),size(xx)); u2 = reshape(ug(Mt+(1:Mt)),size(xx));
p = reshape(pg,size(xx));
if expt=='t'                            % show errors vs known soln...
  i=ceil(ny/2); j=ceil(nx/4); % index interior pt to get pres const
  p = p - p(i,j) + peg(i,j);            % shift const in pres to match known
  eg2 = sum([u1(:)-ue1(:),u2(:)-ue2(:)].^2,2); % squared ptwise vector L2 errs
  if v, figure; subplot(1,2,1); imagesc(gx,gy,log10(reshape(eg2,size(xx)))/2);
    axis xy equal tight; caxis([-16 0]); colorbar; hold on; plot(U.x,'k.-');
    plot(D.x,'k.-'); title('peri vel Stokes BVP: log_{10} u err')
    subplot(1,2,2); imagesc(gx,gy,log10(abs(p-reshape(peg,size(xx)))));
    axis xy equal tight; caxis([-16 0]); colorbar; hold on;
    plot(U.x,'k.-'); plot(D.x,'k.-'); title('log_{10} p err (up to const)'); end
else                                     % just show (velocity,pressure) soln...
  if v, figure; p0 = p(35,1); % pressure const
    contourf(gx,gy,p.*di, p0+pgro*(0:0.05:1)); hold on;
    u0 = 0.1; u1c = min(max(u1,-u0),u0); u2c = min(max(u2,-u0),u0); % clip vel
    colorbar; axis xy equal; quiver(gx,gy, u1c.*di,u2c.*di, 3);
    plot(U.x,'k.-'); plot(D.x,'k.-'); axis([0 uc.d -pi pi]);
    title('peri no-slip p-driven Stokes BVP: (u,p)'), end
end

%keyboard % don't forget to use dbquit to finish otherwise trapped in debug mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function s = quadr(s, N)  % set up periodic trapezoid quadrature on a segment
% Note sign change in normal vs periodicdirpipe.m. 
% Barnett 4/21/13 from testlapSDevalclose.m
t = (1:N)'/N*2*pi; s.x = s.Z(t); s.sp = abs(s.Zp(t)); s.nx = -1i*s.Zp(t)./s.sp;
s.cur = -real(conj(s.Zpp(t)).*s.nx)./s.sp.^2; s.w = 2*pi/N*s.sp; % speed weights
s.t = t; s.cw = 1i*s.nx.*s.w;  % complex weights (incl complex speed)

function h = showseg(s, uc) % plot a segment & possibly its nei copies
if nargin<2, uc.nei = 0; uc.d = 0; end % dummy uc
if uc.nei>0, uc.nei = 1; end % don't plot a ridiculous # of copies
for i=-uc.nei:uc.nei
  plot(s.x+uc.d*i, 'b.-'); axis equal xy;
  l=0.3; hold on; plot([s.x+uc.d*i, s.x+l*s.nx+uc.d*i].', 'k-');
end

% Following kernel matrix fillers are taken from (debugged in) testkernels.m:

function [A T] = SLPmatrix(t,s,mu,a) % single-layer 2D Stokes kernel vel matrix
% Returns 2N-by-2N matrix from src force vector to 2 flow component
% t = target seg (x cols), s = src seg, a = optional translation of src seg
% No option for self-int, gives Inf on diag.   3/2/14
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

function A = SLPpresmatrix(t,s,mu,a) % single-layer 2D Stokes kernel press mat
% Returns N-by-2N matrix from src force vector to pressure value, no self-int.
% t = target seg (x cols), s = src seg, a = optional transl of src seg. 2/1/14
if nargin<4, a = 0; end
N = numel(s.x); M = numel(t.x);
r = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
irr = 1./(conj(r).*r);    % 1/r^2
d1 = real(r); d2 = imag(r);
A = [d1.*irr/(2*pi), d2.*irr/(2*pi)];                  % pressure
A = A .* repmat([s.w(:)' s.w(:)'], [M 1]);           % quadr wei

function [A T] = DLPmatrix(t,s,mu,a) % double-layer 2D Stokes vel kernel matrix
% Returns 2N-by-2N matrix from src force vector to 2 flow components
% If detects self-int, does correct diagonal limit, no jump condition (PV int).
% t = target seg (x cols), s = src seg, a = optional transl of src seg.
% 2nd output is optional traction matrix, without self-eval. 2/1/14, 3/2/14
if nargin<4, a = 0; end
N = numel(s.x); M = numel(t.x);
r = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
irr = 1./(conj(r).*r);    % 1/R^2
d1 = real(r); d2 = imag(r);
rdotny = d1.*repmat(real(s.nx)', [M 1]) + d2.*repmat(imag(s.nx)', [M 1]);
rdotnir4 = rdotny.*(irr.*irr); if nargout<=1, clear rdotny; end
A12 = (1/pi)*d1.*d2.*rdotnir4;  % off diag vel block
A = [(1/pi)*d1.^2.*rdotnir4,   A12;                     % Ladyzhenzkaya
     A12,                      (1/pi)*d2.^2.*rdotnir4];
if numel(s.x)==numel(t.x) & max(abs(s.x+a-t.x))<1e-14
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
  
function A = DLPpresmatrix(t,s,mu,a) % double-layer 2D Stokes kernel press mat
% Returns N-by-2N matrix from src force vector to pressure values, no self-int
% t = target seg (x cols), s = src seg, a = optional transl of src seg. 2/1/14
if nargin<4, a = 0; end
N = numel(s.x); M = numel(t.x);
r = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
irr = 1./(conj(r).*r);    % 1/R^2
d1 = real(r); d2 = imag(r);
rdotn = d1.*repmat(real(s.nx)', [M 1]) + d2.*repmat(imag(s.nx)', [M 1]);
rdotnir4 = rdotn.*(irr.*irr); clear rdotn
A = [(mu/pi)*(-repmat(real(s.nx)', [M 1]).*irr + 2*rdotnir4.*d1), ...
     (mu/pi)*(-repmat(imag(s.nx)', [M 1]).*irr + 2*rdotnir4.*d2)  ];
A = A .* repmat([s.w(:)' s.w(:)'], [M 1]);           % quadr wei

function i = diagind(A,s) %return indices of (shifted) diagonal of square matrix
if nargin<2, s=0; end % if present, s shifts diagonal cyclicly
N = size(A,1); i = sub2ind(size(A), 1:N, mod(s+(0:N-1),N)+1);
