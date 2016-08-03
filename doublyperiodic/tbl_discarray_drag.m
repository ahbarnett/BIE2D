function tbl_discarray_drag
% Generate table of dimensionless drag of regular array of discs, to match
% Table 1 of Greengard-Kropinski, J. Engs. Math. 48: 157–170, 2004.
% Single inclusion (K=1), native quad for matrix fill, ELS.
% Adapted from fig_stoconvK1.m
% Barnett 8/2/16

warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square, more stable

mu = 0.8;   % overall viscosity const for Stokes PDE, can be anything positive
sd = [1 1];  % layer potential representation prefactors for SLP & DLP resp.
jumps = [1 0]; % pressure jumps across R-L and T-B

U.e1 = 1; U.e2 = 1i;     % unit cell; nei=1 for 3x3 scheme
U.nei = 1; [tx ty] = meshgrid(-U.nei:U.nei); U.trlist = tx(:)+1i*ty(:);
m = 22; [U L R B T] = doublywalls(U,m);
proxyrep = @StoSLP;      % sets proxy pt type via a kernel function call
Rp = 1.4; M = 80; % 80;       % proxy params
p.x = Rp * exp(1i*(0:M-1)'/M*2*pi); p = setupquad(p);  % proxy pts

cs = [0.05 0.1:0.1:0.7 0.75 0.76 0.77];  % vol fracs, from G-K 04 paper
Ns = [60*ones(1,5) 100 150 350 800 1200 2000];
close = 0;   % if 0, plain Nystrom for Aelse; if 1, close-eval (slow 14s N=200!)
%cs = 0.75; Ns = 150; close = 1; s.a = 0;   % testing, v slow & only 1e-9 acc!
cs = kron(cs,[1 1]); Ns = kron(Ns,[1 1]) + kron(1+0*Ns,[0 50]); % check all conv
ts = 0*cs; Dcalcs = 0*cs;               % what we compute
for i=1:numel(cs)            % -------------------- loop over filling fractions
  N = Ns(i); r0 = sqrt(cs(i)/pi); h=2*pi*r0/N; delta=(1-2*r0)/h;  % dist
  s = wobblycurve(r0,0,1,N);
  % obstacle no-slip & pressure-drop driving...
  rhs = [zeros(2*N,1); zeros(2*m,1);jumps(1)+0*L.x;zeros(4*m,1);jumps(2)+0*B.x];
  tic
  E = ELSmatrix(s,p,proxyrep,mu,sd,U,close);          % fill (w/ close option)
  co = linsolve(E,rhs,lso);                           % direct bkw stable solve
  ts(i) = toc;
  %fprintf('resid norm = %.3g\n',norm(E*co - rhs))
  sig = co(1:2*N); psi = co(N+1:end);
  %fprintf('density norm = %.3g, proxy norm = %.3g\n',norm(sig), norm(psi))
  J = evalfluxes(s,p,proxyrep,mu,sd,U,co);
  %fprintf('fluxes (%.16g,%.16g)\n',J(1),J(2))
  Dcalcs(i) = 1/(J(1)*mu);                % dimless drag;  note force = 1
  fprintf('c=%g\tr=%.3g\tN=%d (%.2gh)\t%.3gs\tD = %.14e\n',cs(i),r0,N,delta,ts(i),Dcalcs(i))
end                         % ---------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [E A B C Q] = ELSmatrix(s,p,proxyrep,mu,sd,U,close)
% builds matrix blocks for Stokes extended linear system, S+D rep w/ Kress self
% & with close option for A block 3x3 neighbors. Barnett 8/2/16
N = numel(s.x);
if ~close       % plain Nystrom for nei interactions...
  A = sd(1)*srcsum(@StoSLP,U.trlist,[],s,s,mu) + sd(2)*(eye(2*N)/2 + srcsum(@StoDLP,U.trlist,[],s,s,mu));   % notes: DLP gives ext JR term; srcsum self is ok
else            % Nystrom self & close-eval for nei interactions...
  A = sd(1)*StoSLP(s,s,mu) + sd(2)*(eye(2*N)/2 + StoDLP(s,s,mu));  % self w/ JR
  notself = U.trlist(U.trlist~=0);
  A = A + sd(1)*srcsum2(@StoSLP_closeglobal,notself,[],s,s,mu,[],'e') + sd(2)*srcsum2(@StoDLP_closeglobal,notself,[],s,s,mu,[],'e');  % dens=[], 'e' exterior
end
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
