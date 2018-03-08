function stoK1_subdrop
% version of fig_stoconvK1, doubly-periodic 2D pressure-driven Stokes solver.
% Single inclusion (K=1), native quad for matrix fill, close eval for soln.
% Try subtraction of pressure drop via "scattering" split to u = v + w,
% with w solved by GMRES on a mocked-up perifmm.
%
% This version works by using Wen's trick of adding 2 unknown constant
% translational flows to the perifmm representation; this accounts for
% whatever constant offset convention the perifmm has. We also fixed a
% convention for the perifmm which handles inconsistent total force and
% flow by distributing it uniformly (via Z) onto the discrepancy. This
% has no instability. Variation in the solution of the Q solve (empty BVP)
% is removed by enforcing orthogonality R'*xi = 0_3 (if this is not done,
% gmres only gets 6 digits due to xi hence offsets flopping around).
%
% Barnett 1/3/18-1/12/18. 3/6/18 again trying find perifmm form - works.
%
% Todo: rewrite the final eval/flux to use the mock perifmm too.

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
rhs = [zeros(2*N,1); PC]; % driving for ELS
[E,A,B,C,Q] = ELSmatrix(s,p,proxyrep,mu,sd,U);                % fill blocks
%S = svd(E); disp('last few sing vals of E:'); disp(S(end-5:end))
%[~,S,V] = svd(E); figure; plot(V(:,end)); % show the 1 vec in Nul E: all xi
% figure; subplot(2,1,1); imagesc(V); subplot(2,1,2); semilogy(diag(S)); % all low sing vals have r sing vecs assoc w/ xi
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

% build all stacks of 3 "ones-vectors"...
% Y is 2N*3, enforces compat (includes quadr wei). Signs so W'*C=Y' ...
Y = [-s.w' zeros(1,N); zeros(1,N) -s.w'; s.w'.*real(s.nx)' s.w'.*imag(s.nx)']';
Y2 = Y(:,1:2);   % just the force (not incompress) compat
% H is 2N*3, densities generating incompat
H = [ones(1,N) zeros(1,N);zeros(1,N) ones(1,N);real(s.nx)' imag(s.nx)']';
% R is 2M*3, is proxy coeffs proj onto Nul Q
R = [[ones(M,1);zeros(M,1)],[zeros(M,1);ones(M,1)],[real(p.x);imag(p.x)]]/M;
% Z, W are 8m*3;  Z is const discrep incons vecs, W lin func defining compat
w = U.L.w'; W = [[zeros(2*m,1);w;zeros(3*m,1);w;zeros(m,1)],...
                 [zeros(3*m,1);w;zeros(3*m,1);w],...
                 [w;zeros(4*m,1);w;zeros(2*m,1)]]; Z = (W>0); % W'*Z = 2I_3

% solve for v via correction to const sigma: Stokes in U\Omega, has inhomog PC
perim = sum(s.w);
constsig = -kron(jumps',ones(N,1)) / perim / sd(1);  % gives left force = jumps
PC2 = PC - C*constsig;        % PC now consistent
xi = linsolve(Q,PC2,lso);
fprintf('v solve: |xi|=%.3g, resid nrm = %.3g\n',norm(xi), norm(Q*xi-PC2))
fprintf('check v rep matches PC: %.3g\n',norm(C*constsig + Q*xi - PC))
vdata = A*constsig + B*xi;   % v on bdry (exterior limit), rep in 2 basis types

% solve for w
wBC = -vdata;       % RHS for the density-only lin sys
%fprintf('|wBC|_2 = %.3g\nY''*vdata:\n',norm(wBC)), disp(Y'*vdata)
% note v zero-flux but doesn't have zero total vector
%co = randn(2*N,1); [y,wxi] = mockperifmm(co,A,B,C,Q,s,p,U); stop % test mock fmm
%for i=1:3, norm(mockperifmm(H(:,i),A,B,C,Q,s,p,U)), end, stop  % check zero

fprintf('\tis W^TC = Y^T?  yes to %.3g\n',norm(W'*C - Y'))
% test perifmm nullspace is H - it is:
%for i=1:3, norm(mockperifmm(H(:,i),A,B,C,Q,s,p,U)), end
%inv(Y'*H)   % exists! so Y has full rank proj onto H=Nul(Aper)
%fprintf('is b in ran(Aper)?');
fprintf('\t norm(QR) = %.3g\n',norm(Q*R))


[Aper X] = fillAper(A,B,C,Q,s,p,U); % get full Aper matrix
S = svd(Aper); disp('last few sing vals of Aper:'); disp(S(end-5:end)) %cond(Aper)
T = [ones(N,1) zeros(N,1);zeros(N,1) ones(N,1)];  % 2N*2 stack of const flows
al = [0;0];  % default transl flow coeffs, will get solved for too.
if 0       % direct, new, is just wrong
  wsig = linsolve(Aper,wBC,lso);       % direct solve
  fprintf('direct: ||wsig||=%.3g, resid=%.3g\n',norm(wsig),norm(Aper*wsig-wBC))
  wxi = X*wsig;  % create xi aux rep for this density soln
elseif 0       % direct 2-expanded (Wen style) ... works
  Apere = [Aper, T; Y2', zeros(2,2)];      %cond(Apere)
  co = linsolve(Apere,[wBC;0;0],lso);       % direct solve
  wsig = co(1:2*N); al = co(end-1:end);    % get sigma and transl flow coeffs
  fprintf('direct: ||wsig||=%.3g, resid=%.3g\n',norm(wsig),norm(Aper*wsig+T*al-wBC))
  wxi = X*wsig;  % create xi aux rep for this density soln
elseif 1   % ...or, iter for 2-expanded case (Wen style) ... works
  matvec = @(x) mockperifmm(x,A,B,C,Q,s,p,U);   % the perifmm
  matvece = @(x) [matvec(x(1:2*N))+T*x(end-1:end); Y2'*x(1:2*N)]; %expand matvec
  [co,flag,relres,iter,resvec] = gmres(matvece,[wBC;0;0],[],1e-14,2*N);
  its=iter(2),resvec
  % since lost it, re-grab wxi the aux coeffs for soln dens wsig...
  wsig = co(1:2*N); al = co(end-1:end);    % get sigma and transl flow coeffs
  [~,wxi] = mockperifmm(wsig,A,B,C,Q,s,p,U);
else       % ...or, ELS solve for just w part...  works, of course
  wrhs = [wBC;zeros(8*m,1)]; co = linsolve(E,wrhs,lso);  % rhs = [-vdata;0]
  fprintf('for w: ELS resid nrm = %.3g\n',norm(E*co - wrhs))
  wsig = co(1:2*N); wxi = co(2*N+1:end);
end
disp('Y''*wsig (should all be zero):'); disp(Y'*wsig)
%disp('H''*wsig :'); disp(H'*wsig)
%disp('WT.C.wsig'); disp(W'*C*wsig);
%wxi = -Q\(C*wsig);
%disp('||Q.wxi+C.wsig||'); norm(Q*wxi+C*wsig)

disp('R''*wxi :'); disp(R'*wxi)
fprintf('w-solve 1st row chk: %.3g\n',norm(A*wsig + B*wxi + T*al - wBC))
fprintf('w-solve 2nd row chk: %.3g\n',norm(C*wsig + Q*wxi))
P = Z*inv(W'*Z)*W'; Cproj = C-P*C; fprintf('w-solve 2nd row Cproj chk: %.3g\n',norm(Cproj*wsig + Q*wxi))
%norm(P*C*wsig)  % zero when wsig correct

sig = wsig + constsig; psi = wxi + xi;   % combine u = w + v  via their coeffs
co = [sig;psi];                       % unfolded coeffs for total eval

fprintf('density norm = %.3g, proxy norm = %.3g\n',norm(sig), norm(psi))
fprintf('body force + jumps = (%.3g,%.3g)  should vanish\n',s.w'*sig(1:N)+abs(U.e2)*jumps(1),s.w'*sig(N+1:end)+abs(U.e1)*jumps(2))
% todo: replace this eval and fluxes w/ mock perifmm call...
[u p0] = evalsol(s,p,proxyrep,mu,sd,U,z,co);        % native quad (far) test
u = u + al;         % note: the rep for soln u includes transl constants
fprintf('u at pt = (%.16g,%.16g)  \t(est abs err: %.3g)\n',u(1),u(2),norm(u-uex))
J = evalfluxes(s,p,proxyrep,mu,sd,U,co);
J = J + al;
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

function [y xi] = mockperifmm(x,A,B,C,Q,s,p,U)  % apply Aper for now
% x is inclusion density vec. Rest are the dense ELS blocks, for whatever rep.
% homog periodicity.
% Outputs:
% y = Aper*x, as Stokes 2d2p perifmm would do. xi is the aux coeffs used.
% Barnett 1/10/18; retry 3/7/18
%[Aper genxi] = fillAper(A,B,C,Q,s,p); y=Aper*x; xi=genxi*x; return % test!
M = numel(p.x); N = numel(s.x); m=numel(U.L.x);   % # proxies, bdry, wall nodes
% build the "ones-vectors"...
% Y is 2N*3, enforces compat (includes quadr wei)
Y = [s.w' zeros(1,N); zeros(1,N) s.w'; s.w'.*real(s.nx)' s.w'.*imag(s.nx)']';
% H is 2N*3, densities generating incompat
H = [ones(1,N) zeros(1,N);zeros(1,N) ones(1,N);real(s.nx)' imag(s.nx)']';
% R is 2M*3, is proxy
R = [[ones(M,1);zeros(M,1)],[zeros(M,1);ones(M,1)],[real(p.x);imag(p.x)]]/M;
% Z, W are 8m*3;  Z is const discrep incons vecs, W lin func defining compat
w = U.L.w'; W = [[zeros(2*m,1);w;zeros(3*m,1);w;zeros(m,1)],...
                 [zeros(3*m,1);w;zeros(3*m,1);w],...
                 [w;zeros(4*m,1);w;zeros(2*m,1)]]; Z = (W>0); % W'*Z = 2I_3
P = Z*inv(W'*Z)*W';
%rng(0); R = randn(size(R));  % try some other fixed R: only gets 7 digits! why?
Qtilde = Q + Z*R';   % since dim nul Q = 3, enforces a convention R'*xi = 0_3
% Notice how GMRES only gets 6 digits if such a convention is not applied
% (xi flops around in the 3d-nullspace and messes up convergence)
discrep = C*x;
g = - (discrep - P*discrep);  % project to a consistent discrep for empty BVP
xi = Qtilde \ g;    % solve aux dofs;  checked R'*xi is indeed 0_3
y = A*x + B*xi;   % add two contribs for rep
% todo: also output J1,J2 since will want in real perifmm.

function [Aper X] = fillAper(A,B,C,Q,s,p,U)  % densely fill Aper
% A..Q are the dense ELS blocks, for whatever rep.
% Outputs: Aper, as Stokes 2d2p perifmm would apply. if til, perturb Q,B->tilde.
% genxi is mat which gives xi from the Aper soln, ie -X in DPLS.
% Barnett 3/7/18
M = numel(p.x); N = numel(s.x); m=numel(U.L.x);  % # proxies, bdry & wall nodes
% build the "ones-vectors"...
% H is 2N*3, densities generating incompat
H = [ones(1,N) zeros(1,N);zeros(1,N) ones(1,N);real(s.nx)' imag(s.nx)']';
% R is 2M*3, is proxy
R = [[ones(M,1);zeros(M,1)],[zeros(M,1);ones(M,1)],[real(p.x);imag(p.x)]]/M;
% Z, W are 8m*3;  Z is const discrep incons vecs, W lin func defining compat
w = U.L.w'; W = [[zeros(2*m,1);w;zeros(3*m,1);w;zeros(m,1)],...
                 [zeros(3*m,1);w;zeros(3*m,1);w],...
                 [w;zeros(4*m,1);w;zeros(2*m,1)]]; Z = (W>0); % W'*Z = 2I_3
P = Z*inv(W'*Z)*W';
%R = 0*R; R(1:3,1:3)=eye(3);  % arb fixed R
%rng(0); R = randn(size(R));  % try arb but fixed R - works equally well here
Qtilde = Q + Z*R';   % enforce a convention R'*xi = 0
G = -C + P*C;
X = Qtilde \ G;    % matrix version of xi
fprintf('\tnorm X = %.3g\t\tnorm(R^TX) = %.3g\n',norm(X),norm(R'*X))
%RTX = R'*X; RTX(:,end/2)  % if don't use P, then only 2nd cmpt (row) = 0.
%inv(W'*Z) %  I/2
Aper = A + B*X;
