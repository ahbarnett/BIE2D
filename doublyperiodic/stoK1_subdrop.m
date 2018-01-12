function stoK1_subdrop
% version of fig_stoconvK1, doubly-periodic 2D pressure-driven Stokes solver.
% Single inclusion (K=1), native quad for matrix fill, close eval for soln.
% Try subtraction of pressure drop via "scattering" split to u = v + w,
% with w solved by GMRES on mocked-up perifmm.
% Barnett 1/3/18-1/12/18

disp('2D 2-periodic Stokes pressure-driven, K=1 inclusion, via split method')
warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square

mu = 0.7;   % overall viscosity const for Stokes PDE
sd = [1 1];  % layer potential representation prefactors for SLP & DLP resp.
jumps = [1 0];  % pressure drops as go L-to-R and T-to-B

U.e1 = 1; U.e2 = 1i;     % unit cell; nei=1 for 3x3 scheme
U.nei = 1; [tx ty] = meshgrid(-U.nei:U.nei); U.trlist = tx(:)+1i*ty(:);
m = 22; U = doublywalls(U,m);
proxyrep = @StoSLP;      % sets proxy pt type via a kernel function call
Rp = 1.4; M = 70;  % proxy params - 80 is 1-2 digits worse than 70 or 120 - why?
p.x = Rp * exp(1i*(0:M-1)'/M*2*pi); p = setupquad(p);  % proxy pts

a = 0.7; b = 0.15;   % worm params, spills horizontally out of any unit cell
z = .1+.4i;                                         % a test pt
uex = [0.016778793238;0.005152952237]; flux1ex = 0.008234042360; %for a=.7,b=.15

N = 150; s = wormcurve(a,b,N); s.a = mean(s.x);  % a needed for p ext close
%s.x = s.x + 0.1i;  % check translational invariance of flux

PC = [zeros(2*m,1);jumps(1)+0*U.L.x;zeros(4*m,1);jumps(2)+0*U.B.x]; % physical
rhs = [zeros(2*N,1); PC]; % driving
[E,A,B,C,Q] = ELSmatrix(s,p,proxyrep,mu,sd,U);                % fill blocks
%figure; imagesc(E); colorbar; drawnow
if 0   % direct ELS solve (no scatt split), checks the setup is good.
  co = linsolve(E,rhs,lso);
  fprintf('ELS resid nrm = %.3g\n',norm(E*co - rhs))
  sig = co(1:2*N); psi = co(2*N+1:end);
  fprintf('body force + jumps = (%.3g,%.3g)  should vanish\n',s.w'*sig(1:N)+abs(U.e1)*jumps(1),s.w'*sig(N+1:end)+abs(U.e2)*jumps(2))
  [u p0] = evalsol(s,p,proxyrep,mu,sd,U,z,co);        % native quad (far) test
  fprintf('u at pt = (%.16g,%.16g)  \t(est abs err: %.3g)\n',u(1),u(2),norm(u-uex))
  J = evalfluxes(s,p,proxyrep,mu,sd,U,co);
  fprintf('fluxes (%.16g,%.16g)   (est abs err: %.3g)\n',J(1),J(2),abs(J(1)-flux1ex))
end
  
% solve for v via correction to const sigma: Stokes in U\Omega, has inhomog PC
perim = sum(s.w);
constsig = -kron(jumps',ones(N,1)) / perim / sd(1);  % gives left force = jumps
PC = PC - C*constsig;        % PC now consistent
xi = linsolve(Q,PC,lso);
fprintf('v solve: |xi|=%.3g, resid nrm = %.3g\n',norm(xi), norm(Q*xi-PC))
vdata = A*constsig + B*xi;   % v on bdry (exterior limit), rep in 2 basis types

% build the "ones-vectors"...
% Y is 2N*3, enforces compat (includes quadr wei)
Y = [s.w' zeros(1,N); zeros(1,N) s.w'; s.w'.*real(s.nx)' s.w'.*imag(s.nx)']';
% H is 2N*3, densities generating incompat
H = [ones(1,N) zeros(1,N);zeros(1,N) ones(1,N);real(s.nx)' imag(s.nx)']';
% R is 2M*3, is proxy coeffs proj onto Nul Q
R = [[ones(M,1);zeros(M,1)],[zeros(M,1);ones(M,1)],[real(p.x);imag(p.x)]]/M;

% solve for w
wBC = -vdata;       % RHS for the density-only lin sys
%fprintf('|wBC|_2 = %.3g\nY''*vdata:\n',norm(wBC)), disp(Y'*vdata)
% note v zero-flux but doesn't have zero total vector
%co = randn(2*N,1); [y,wxi] = mockperifmm(co,A,B,C,Q,s,p); stop % test mock fmm
%for i=1:3, norm(mockperifmm(H(:,i),A,B,C,Q,s,p)), end, stop  % check zero

til = 1;      % 0: don't perturb B,Q; 1 use Btilde,Qtilde pert via HR'
[Aper genxi] = fillAper(A,B,C,Q,s,p,til); % get full Aper matrix
if 0  % direct
  S = svd(Aper); disp('last few sing vals of Aper:'); disp(S(end-5:end))
  wsig = linsolve(Aper,wBC,lso);       % direct solve
  fprintf('direct: ||wsig||=%.3g, resid=%.3g\n',norm(wsig),norm(Aper*wsig-wBC))
  wxi = genxi*wsig;  % create xi aux rep for this density soln
  if til, wsig = wsig + H*(R'*wxi); end % correct the density
  
  %Y = [s.w' zeros(1,N); zeros(1,N) s.w'; s.w'.*real(s.nx)' s.w'.*imag(s.nx)']';
  %H = [ones(1,N) zeros(1,N);zeros(1,N) ones(1,N);real(s.nx)' imag(s.nx)']';
  %al = -(Y'*H)\(Y'*wsig)
  %wsig = wsig + H*al; % project out wsig soln
elseif 1  % ...or, iter
  matvec = @(x) mockperifmm(x,A,B,C,Q,s,p,til);
  [wsig,flag,relres,iter,resvec] = gmres(matvec,wBC,[],1e-14,2*N);
  %[wsig,flag,relres,iter,resvec] = gmres(matvec,wBC,[],1e-14,2*N,[],[],randn(2*N,1));  % random start vec x_0
  its=iter(2),resvec
  % since lost it, re-grab wxi the aux coeffs for soln dens wsig...
  [~,wxi] = mockperifmm(wsig,A,B,C,Q,s,p,til);
  if til, wsig = wsig + H*(R'*wxi); end % correct the density
else    % ...or, ELS solve for just w part...  works, of course
  wrhs = [wBC;zeros(8*m,1)]; co = linsolve(E,wrhs,lso);  % rhs = [-vdata;0]
  fprintf('for w: ELS resid nrm = %.3g\n',norm(E*co - wrhs))
  wsig = co(1:2*N); wxi = co(2*N+1:end);
end
disp('Y''*wsig (should all be zero):'); disp(Y'*wsig)

sig = wsig + constsig; psi = wxi + xi;   % combine u = w + v  via their coeffs
co = [sig;psi];                       % unfolded coeffs for total eval

fprintf('density norm = %.3g, proxy norm = %.3g\n',norm(sig), norm(psi))
fprintf('body force + jumps = (%.3g,%.3g)  should vanish\n',s.w'*sig(1:N)+abs(U.e2)*jumps(1),s.w'*sig(N+1:end)+abs(U.e1)*jumps(2))
[u p0] = evalsol(s,p,proxyrep,mu,sd,U,z,co);        % native quad (far) test
fprintf('u at pt = (%.16g,%.16g)  \t(est abs err: %.3g)\n',u(1),u(2),norm(u-uex))
J = evalfluxes(s,p,proxyrep,mu,sd,U,co);
fprintf('fluxes (%.16g,%.16g)   (est abs err: %.3g)\n',J(1),J(2),abs(J(1)-flux1ex))

%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u p] = evalsol(s,pr,proxyrep,mu,sd,U,z,co,close) % eval soln rep u,p
% z = list of targets as C values. pr = proxy struct, s = source curve struct
% co = full coeff vec, U = unit cell struct, mu = viscosity
% sd = prefactors for source rep, close enables special close quadr.
% Note: not for self-eval since srcsum2 used, 6/30/16
if nargin<9, close=0; end              % default is plain native quadr
z = struct('x',z);                     % make targets z a segment struct
N = numel(s.x);
sig = co(1:2*N); psi = co(2*N+1:end);  % split into sig (density) & psi (proxy)
if close, S = @(t,s,mu,dens) StoSLP_closeglobal(t,s,mu,dens,'e'); % exterior
  D = @(t,s,mu,dens) StoDLP_closeglobal(t,s,mu,dens,'e');
else, S = @StoSLP; D = @StoDLP; end    % NB now native & close have same 4 args
if nargout==1                          % don't want pressure output
  u = proxyrep(z,pr,mu,psi);           % init sol w/ proxies (always far)
  u = u + sd(1)*srcsum2(S,U.trlist,[],z,s,mu,sig) + sd(2)*srcsum2(D,U.trlist,[],z,s,mu,sig);
else  
  [u p] = proxyrep(z,pr,mu,psi);       % init sol w/ proxies (always far)
  [uS pS] = srcsum2(S,U.trlist,[],z,s,mu,sig);
  [uD pD] = srcsum2(D,U.trlist,[],z,s,mu,sig);
  u = u + sd(1)*uS + sd(2)*uD;
  p = p + sd(1)*pS + sd(2)*pD;
end

function J = evalfluxes(s,p,proxyrep,mu,sd,U,co)
% inputs as in evalsol. Uses Gary-inspired bdry of 3x3 block far-field method
if U.nei~=1, warning('U.neu must equal 1'); end
w = U.L.w; if norm(w-U.B.w)>1e-14, error('L and B must have same weights'); end
m = numel(w);
N = numel(s.x); sig = co(1:2*N); psi = co(2*N+1:end);
t.x = [U.L.x;U.B.x]; t.nx = [U.L.nx;U.B.nx];  % 2-wall target
v = proxyrep(t, p, mu, psi);                  % flow vel on L+B
u = v(1:2*m).*real(t.nx) + v(2*m+1:end).*imag(t.nx); % proxy contrib to u = v.n
J = sum(reshape(u,[m 2])'.*([1;1]*w),2);         % .... do its quadr on L,B 
for i=-1:1      % set up big loop of 12 walls, just their nodes
  x{i+2} = U.L.x - U.e1 +i*U.e2; x{i+5} = x{i+2} + 3*U.e1;
  x{i+8} = U.B.x +i*U.e1 -U.e2; x{i+11} = x{i+8} + 3*U.e2;
end
t.x = vertcat(x{:}); t.nx = [repmat(U.L.nx,[6 1]); repmat(U.B.nx,[6 1])];
v = sd(1)*StoSLP(t,s,mu,sig) + sd(2)*StoDLP(t,s,mu,sig);  % ctr copy only
u = v(1:12*m).*real(t.nx) + v(12*m+1:end).*imag(t.nx);   % v.n, as col vec
amts = [0 0 0 3 3 3 -1 -2 -3 1 2 3; -1 -2 -3 1 2 3 0 0 0 3 3 3];  % wall wgts
J = J + sum(repmat(u',[2 1]).*kron(amts,w),2);   % weight each wall

function [E A B C Q] = ELSmatrix(s,p,proxyrep,mu,sd,U)
% builds matrix blocks for Stokes extended linear system, S+D rep w/ Kress self
N = numel(s.x);
A = sd(1)*srcsum(@StoSLP,U.trlist,[],s,s,mu) + sd(2)*(eye(2*N)/2 + srcsum(@StoDLP,U.trlist,[],s,s,mu));   % notes: DLP gives exterior JR term; srcsum self is ok
B = proxyrep(s,p,mu);     % map from proxy density to vel on curve
C = Cblock(s,U,mu,sd);
[QL,~,QLt] = proxyrep(U.L,p,mu); [QR,~,QRt] = proxyrep(U.R,p,mu); % vel, tract
[QB,~,QBt] = proxyrep(U.B,p,mu); [QT,~,QTt] = proxyrep(U.T,p,mu);
Q = [QR-QL; QRt-QLt; QT-QB; QTt-QBt];
E = [A B; C Q];

function C = Cblock(s,U,mu,sd)     % fill C from source curve s to U walls
% sd controls prefactors on SLP & DLP for Stokes rep on the obstacle curve
n = U.nei; e1 = U.e1; e2 = U.e2; S = @StoSLP; D = @StoDLP;      % abbrevs
trlist = n*e1 + (-n:n)*e2;
[CLS,~,TLS] = srcsum(S,trlist,[],U.L,s,mu);
[CLD,~,TLD] = srcsum(D,trlist,[],U.L,s,mu);
trlist = -n*e1 + (-n:n)*e2;
[CRS,~,TRS] = srcsum(S,trlist,[],U.R,s,mu);
[CRD,~,TRD] = srcsum(D,trlist,[],U.R,s,mu);
trlist = (-n:n)*e1 + n*e2;
[CBS,~,TBS] = srcsum(S,trlist,[],U.B,s,mu);
[CBD,~,TBD] = srcsum(D,trlist,[],U.B,s,mu);
trlist = (-n:n)*e1 - n*e2;
[CTS,~,TTS] = srcsum(S,trlist,[],U.T,s,mu);
[CTD,~,TTD] = srcsum(D,trlist,[],U.T,s,mu);
C = sd(1)*[CRS-CLS; TRS-TLS; CTS-CBS; TTS-TBS]  +...
    sd(2)*[CRD-CLD; TRD-TLD; CTD-CBD; TTD-TBD];

function [zz ii] = extgrid(gx,gy,s,U)  % grid points and indices outside Omegas
% given gx,gy 1d x,y grids, s=curve segment, U = unit cell struct, return
% zz = col list of grid pts as C-#s, and ii = logical index array if outside
% Omega and all images
[xx yy] = meshgrid(gx,gy); zz = (xx(:)+1i*yy(:)); clear xx yy;
ii = true(size(zz));               % indices outside Omega or its copies
for i=-1:1, for j=-1:1, si = s.inside(zz+U.e1*i+U.e2*j); ii(si)=false; end, end

function [y xi] = mockperifmm(x,A,B,C,Q,s,p,til)  % densely apply Aper for now
% x is inclusion density vec. Rest are the dense ELS blocks, for whatever rep.
% First orthog projects to make density compatible, per defn is homog PCs.
% if til=0, leave Q,B alone (default); if til=1 replace by Qtilde,Btilde.
% Outputs:
% y = Aper*x, as Stokes 2d2p perifmm would do. xi is the aux coeffs used.
% Barnett 1/10/18
%[Aper genxi] = fillAper(A,B,C,Q,s,p,til); y=Aper*x; xi=genxi*x; return % test!
if nargin<8, til=0; end
M = numel(p.x); N = numel(s.x);   % # proxies, bdry nodes
% build the "ones-vectors"...
% Y is 2N*3, enforces compat (includes quadr wei)
Y = [s.w' zeros(1,N); zeros(1,N) s.w'; s.w'.*real(s.nx)' s.w'.*imag(s.nx)']';
% H is 2N*3, densities generating incompat
H = [ones(1,N) zeros(1,N);zeros(1,N) ones(1,N);real(s.nx)' imag(s.nx)']';
% R is 2M*3, is proxy
R = [[ones(M,1);zeros(M,1)],[zeros(M,1);ones(M,1)],[real(p.x);imag(p.x)]]/M;

xt=x; %al = -(Y'*H)\(Y'*x); xt = x + H*al;    % project x to xt, so Y'*xt=0_3
g = -C*xt;  % eval the discrep to kill
%fprintf('mockperifmm: |g|_2 = %.3g\n',norm(g))
Qtilde = Q; if til, Qtilde = Q + (C*H)*R'; end
xi = Qtilde\g;   % find aux coeffs to match it. %bare Q amplifies incompat err by 1e9!

%al = -(R'*R)\(R'*xi); xi = xi + R*al;    % project xi so R'*xt=0_3

%fprintf('mockperifmm: |xt|_2 = %.3g,  |xi|_2 = %.3g,  |B.xi|=%.3g\n',norm(xt),norm(xi), norm(B*xi))
Btilde=B; if til, Btilde = B + (A*H)*R'; end
y = A*xt + Btilde*xi;   % "3x3 nr FMM + aux eval". note xi=-"X" in matrix case
%fprintf('mockperifmm: |y|_2 = %.3g\n',norm(y))

%al = -(Y'*Y)\(Y'*y); y = y + Y*al;    % project y so Y'*y=0_3
%Y'*y

function [Aper genxi] = fillAper(A,B,C,Q,s,p,til)  % densely fill Aper for now
% A..Q are the dense ELS blocks, for whatever rep.
% Outputs: Aper, as Stokes 2d2p perifmm would apply. if til, perturb Q,B->tilde.
% genxi is mat which gives xi from the Aper soln, ie -X in DPLS.
% Barnett 1/12/18
M = numel(p.x); N = numel(s.x);   % # proxies, bdry nodes
% build the "ones-vectors"...
% H is 2N*3, densities generating incompat
H = [ones(1,N) zeros(1,N);zeros(1,N) ones(1,N);real(s.nx)' imag(s.nx)']';
% R is 2M*3, is proxy
R = [[ones(M,1);zeros(M,1)],[zeros(M,1);ones(M,1)],[real(p.x);imag(p.x)]]/M;
Qtilde = Q;
Btilde = B; 
if til, Qtilde = Q + (C*H)*R'; Btilde = B + (A*H)*R'; end
X = Qtilde\C;
Aper = A-Btilde*X;
genxi=-X;
