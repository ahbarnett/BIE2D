function discarray_effcond
% Effective conductivity (kappa) of infinite disc array.
% Laplace BVP, single inclusion (K=1), native quad for matrix fill, ELS.
% Adapted from fig_discarray_drag. Uses helsingeffcond.m
% Barnett 9/27/17. Fixed close and reparam 9/28/17

warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square, more stable

jumps = [1 0];   % pressure driving (jumps) across R-L and T-B
schur = 1;       % 0 for rect ELS direct solve; 1 for Schur iterative solve
verbhel = 0;     % print helsing ref soln info

U.e1 = 1; U.e2 = 1i;     % unit cell; nei=1 for 3x3 scheme
U.nei = 1; [tx ty] = meshgrid(-U.nei:U.nei); U.trlist = tx(:)+1i*ty(:);
m = 22; [U L R B T] = doublywalls(U,m);
proxyrep = @LapSLP;      % sets proxy pt type via a kernel function call
Rp = 1.4; M = 80;       % proxy params
p.x = Rp * exp(1i*(0:M-1)'/M*2*pi); p = setupquad(p);  % proxy pts

%cs = 2; sigs = 10; Ns = 120; close=0;  % geom and conductivity params
%cs = 10; sigs = 10; Ns = 240; close=1;  % geom and conductivity params
%cs = 5; sigs = 10; Ns = 900; close=0;  % geom and conductivity params
%cs = 10; sigs = 10; Ns = 800; reparam=1; close=0;  % geom, conductivity params
cs = 100; sigs = 100; Ns = 200; reparam=1; close=1;  % geom, conductivity params
%cs = 1000; sigs = 1000; Ns = 240; reparam=1; close=1;  % geom, conductivity params
%cs = 10; sigs = 10; Ns = 200; close=1;  % geom and conductivity params, fails
% reparam: 0: plain trap rule; 1 reparam bunching near touching pts
% close: if 0, plain Nystrom for Aelse; if 1, close-eval (slow)

chkconv = 0;    % if true, check each ans vs bumped-up N
if chkconv, cs=kron(cs,[1 1]); sigs=kron(sigs,[1 1]); Ns = kron(Ns,[1 1]) + kron(1+0*Ns,[0 400]); end
ts = 0*cs; effconds = 0*cs;               % what we compute
for i=1:numel(cs)            % -------------------- loop over filling fractions
  N = Ns(i); r0 = 0.5*sqrt(1-1/cs(i)^2);   % disc radius
  h=2*pi*r0/N; delta=(1-2*r0)/h;  % quadr spacing h; h-scaled closeness measure
  lam = (sigs(i)-1)/(sigs(i)+1);   % lambda contrast param
  s = wobblycurve(r0,0,1,N); s.a = 0;  % interior pt needed for close only
  if reparam
    be = .7*log(cs(i))/log(2); s = reparam_bunched(s,be); s.cur = (1/r0)*ones(N,1);
  end, %figure; showsegment(s);  % check bunching
  g = [jumps(1)+0*L.x; 0*L.x; jumps(2)+0*B.x; 0*B.x]; rhs = [0*s.x; g]; %driving
  tic
  [E,A,Bm,C,Q] = ELSmatrix(s,p,proxyrep,U,lam,close); % fill (w/ close option)
  if ~schur
    co = linsolve(E,rhs,lso);     % ..... direct bkw stable solve of ELS
    iter=[0 0]; resvec = 0;     % dummies
    tau = co(1:N); psi = co(N+1:end);
  else                               % ..... schur, square well-cond solve
    R = ones(M,1)/M; H = ones(N,1);  % 1s vectors
    Qtilde = Q + (C*H)*R';               % Gary version of Schur
    X = linsolve(Qtilde,C,lso); y = linsolve(Qtilde,g,lso);
    [tau,flag,relres,iter,resvec] = gmres(@(x) A*x - Bm*(X*x), -Bm*y, [], 1e-15, N);  % note iter has 2 elements: 2nd is what want
    xi = y - X*tau;
    co = [tau;xi];  % build full soln vector
  end
  ts(i) = toc;
  %fprintf('resid norm = %.3g\n',norm(E*co - rhs))
  %fprintf('density norm = %.3g, proxy norm = %.3g\n',norm(tau), norm(psi))
  J = evalfluxes(s,p,proxyrep,U,co);
  %fprintf('fluxes (%.16g,%.16g)\n',J(1),J(2))
  effconds(i) = J(1);
  effcondse(i) = helsingeffcond(cs(i),sigs(i),[],verbhel); % reference soln
  fprintf('c=%g r=%.3g\tN=%d (%.2gh)\t%.3gs\t%d its %.3g\tkap=%.14e\n',cs(i),r0,N,delta,ts(i),iter(2),resvec(end),effconds(i))
  fprintf('\t\t\t\trel err from helsing ref = %.3g\n',abs((effconds(i)-effcondse(i))/effcondse(i)))
  if chkconv & mod(i,2)==0, fprintf('\t\t\t\test rel err = %.3g\n',abs((effconds(i)-effconds(i-1))/effconds(i))), end
end                         % ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = evalsol(s,p,proxyrep,U,z,co)
% z = list of targets as C values. p = proxy struct, s = source curve struct
% co = full coeff vec, U = unit cell struct.  Does potential value only.
N = numel(s.x);
tau = co(1:N); psi = co(N+1:end);     % split up solution coeffs (psi=proxy)
u = proxyrep(struct('x',z), p, psi);  % init u w/ proxies eval
for i=1:numel(U.trlist)         % sum potentials faster than srcsum matrices
  u = u + srcsum(@LapSLP, U.trlist(i),[], struct('x',z), s, tau);  % pot of SLP
end

function J = evalfluxes(s,p,proxyrep,U,co)
% inputs as in evalsol. Uses Gary-inspired bdry of 3x3 block far-field method
if U.nei~=1, warning('U.neu must equal 1'); end
w = U.L.w; if norm(w-U.B.w)>1e-14, error('L and B must have same weights'); end
m = numel(w);
N = numel(s.x); tau = co(1:N); psi = co(N+1:end);
t.x = [U.L.x;U.B.x]; t.nx = [U.L.nx;U.B.nx];  % 2-wall target
[~,Tn] = proxyrep(t, p); u = Tn * psi;   % proxy contrib to un only
J = sum(reshape(u,[m 2])'.*([1;1]*w),2);         % .... do its quadr on L,B 
for i=-1:1      % set up big loop of 12 walls, just their nodes
  x{i+2} = U.L.x - U.e1 +i*U.e2; x{i+5} = x{i+2} + 3*U.e1;
  x{i+8} = U.B.x +i*U.e1 -U.e2; x{i+11} = x{i+8} + 3*U.e2;
end
t.x = vertcat(x{:}); t.nx = [repmat(U.L.nx,[6 1]); repmat(U.B.nx,[6 1])];
[~,un] = LapSLP(t, s, tau);          % central density only, many target walls
amts = [0 0 0 3 3 3 -1 -2 -3 1 2 3; -1 -2 -3 1 2 3 0 0 0 3 3 3];  % wall wgts
J = J + sum(repmat(un',[2 1]).*kron(amts,w),2);   % weight each wall

function [E A B C Q] = ELSmatrix(s,p,proxyrep,U,lam,close)
% builds matrix blocks for extended linear system, SLP rep on curve, w/ close
% as option, and lambda conductivity param.
N = numel(s.x)
if ~close          % plain Nystrom for nei interactions...
  [~,A] = srcsum(@LapSLP, U.trlist,[], s,s);   % directly summed self-int matrix
else
  [~,A] = LapSLP(s,s);
  notself = U.trlist(U.trlist~=0);
  [~,Ans] = srcsum2(@LapSLP_closeglobal, notself,[],s,s, [],'e'); % dens=[], 'e' exterior
  A = A + Ans;
end
A = A + eye(N)/(2*lam);        % lambda-modified Neumann exterior jump relation
[~,B] = proxyrep(s,p);         % Neu data from proxies
C = Cblock(s,U,@LapSLP);
[QL QLn] = proxyrep(U.L,p); [QR QRn] = proxyrep(U.R,p);
[QB QBn] = proxyrep(U.B,p); [QT QTn] = proxyrep(U.T,p);
Q = [QR-QL; QRn-QLn; QT-QB; QTn-QBn];
E = [A B; C Q];

function C = Cblock(s,U,densrep)     % fill C from source curve s to U walls
% densrep is handle to LapSLPmat/LapDLPmat depending on density type on curve s
nei = U.nei; N = numel(s.x); m = numel(U.L.x);
[CL CLn] = srcsum(densrep,nei*U.e1 + (-nei:nei)*U.e2,[], U.L,s);
[CR CRn] = srcsum(densrep,-nei*U.e1 + (-nei:nei)*U.e2,[], U.R,s);
[CB CBn] = srcsum(densrep,(-nei:nei)*U.e1 + nei*U.e2,[], U.B,s);
[CT CTn] = srcsum(densrep,(-nei:nei)*U.e1 - nei*U.e2,[], U.T,s);
C = [CR-CL; CRn-CLn; CT-CB; CTn-CBn];
