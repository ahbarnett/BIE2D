function fig_stoconvK1
% make Stokes no-slip periodic convergence and soln plots.
% Single inclusion (K=1), native quad for matrix fill, close eval for soln.
% Adapted from fig_lapconvK1.m
% Barnett 6/7/16. 6/30/16 brought into BIE2D.
% X,y,R,H notation, Gary's V=CH, Alex's nullspace fix, nonrandom. 8/17/16

warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square

mu = 0.7;   % overall viscosity const for Stokes PDE
sd = [1 1];  % layer potential representation prefactors for SLP & DLP resp.
jumps = [1 0]; %[0 1]; % pressure jumps across R-L and T-B (not for known soln)

U.e1 = 1; U.e2 = 1i;     % unit cell; nei=1 for 3x3 scheme
U.nei = 1; [tx ty] = meshgrid(-U.nei:U.nei); U.trlist = tx(:)+1i*ty(:);
m = 22; [U L R B T] = doublywalls(U,m);
proxyrep = @StoSLP;      % sets proxy pt type via a kernel function call
Rp = 1.4; M = 70;  % proxy params - 80 is 1-2 digits worse than 70 or 120 - why?
p.x = Rp * exp(1i*(0:M-1)'/M*2*pi); p = setupquad(p);  % proxy pts

a = 0.7; b = 0.15;   % worm params, spills horizontally out of any unit cell
uex = [0.016778793238;0.005152952237]; flux1ex = 0.008234042360;  % a=.7,b=.15
% known soln stokeslet location & force: must be deep inside Omega (careful)
src.x = 0.42+.23i; src.w = 1; src.nx = 1+0i; fsrc = [0.8;0.6]; src.nei = 5;

% -------------------------- single soln and plot
N = 150; s = wormcurve(a,b,N); s.a = mean(s.x);  % a needed for p ext close
%s.x = s.x + 0.1i;  % check translational invariance of flux
% obstacle no-slip & pressure-drop driving...
rhs = [zeros(2*N,1); zeros(2*m,1);jumps(1)+0*L.x;zeros(4*m,1);jumps(2)+0*B.x];
tic
E = ELSmatrix(s,p,proxyrep,mu,sd,U);                % fill
co = linsolve(E,rhs,lso);                           % direct bkw stable solve
toc
%S = svd(E); disp('last few sing vals of E:'), S(end-5:end) % dim Nul E = 1
fprintf('resid norm = %.3g\n',norm(E*co - rhs))
sig = co(1:2*N); psi = co(2*N+1:end);
fprintf('density norm = %.3g, proxy norm = %.3g\n',norm(sig), norm(psi))
fprintf('body force + jumps = (%.3g,%.3g)  should vanish\n',s.w'*sig(1:N)+abs(U.e1)*jumps(1),s.w'*sig(N+1:end)+abs(U.e2)*jumps(2))
z = .1+.4i;                                         % test pt
[u p0] = evalsol(s,p,proxyrep,mu,sd,U,z,co);        % native quad (far) test
fprintf('u at pt = (%.16g,%.16g)  \t(est abs err: %.3g)\n',u(1),u(2),norm(u-uex))
tic; J = evalfluxes(s,p,proxyrep,mu,sd,U,co); toc
fprintf('fluxes (%.16g,%.16g)   (est abs err: %.3g)\n',J(1),J(2),abs(J(1)-flux1ex))
if 0   % soln flow figure
  nx = 201; gx = 0.5*((0:nx-1)/(nx-1)*2-1); ng = nx^2;  % fine pres eval grid
  gy = gx; [zz ii] = extgrid(gx,gy,s,U);
  pg = nan(ng,1);                           % pressure on fine grid
  tic; [~,pg(ii)] = evalsol(s,p,proxyrep,mu,sd,U,zz(ii),co,1); toc
  pg = reshape(pg,[nx nx]) - p0;            % shift to pres=0 at test pt
  figure; contourf(gx,gy,pg,[-.6:.1:.6]); colormap(jet(256)); hold on;
  nx = 26; gx = 0.5*((0:nx-1)/(nx-1)*2-1); ng = nx^2;  % coarse vel eval grid
  gy = gx; [zz ii] = extgrid(gx,gy,s,U);
  ug = nan(2*ng,1);
  tic; ug([ii;ii]) = evalsol(s,p,proxyrep,mu,sd,U,zz(ii),co,1); toc
  u1 = reshape(ug(1:ng),[nx nx]); u2 = reshape(ug(ng+1:end),[nx nx]); % u cmpts
  % the following sends only the non-nan parts of grid since otherwise get dots:
  quiver(real(zz(ii)),imag(zz(ii)),u1(ii),u2(ii),2.0,'k-');  % show vec field
  n=1; for i=-n:n, for j=-n:n, plot([s.x;s.x(1)]+j+1i*i,'-'); end,end  % curves
  plot([L.x;R.x;T.x;B.x],'b.');  % wall nodes
  axis xy equal off; plot(z,'k.'); axis([0 1 0 1]-0.5);
  %for i=-1:1, for j=-1:1, plot(src.x+U.e1*i+U.e2*j,'k*'); end, end % known
  text(0,0,'$\Omega$','interpreter','latex','fontsize',14);
  text(-.58,.45,'(a)'); %,'fontsize',14);
  drawnow; %set(gcf,'paperposition',[0 0 4 4]); print -depsc2 figs/stosolK1.eps
end

if 1, Ns = 30:10:230;   % ------------------------  N-convergence
us = nan(2,numel(Ns)); ust = us; Js = us; Jst = Js;
es = nan(1,numel(Ns)); res = es; rest = es;
uek = knownsol(U,z,src,fsrc,mu);   % known stokeslet 3x3 grid soln, fixed
%v = zeros(8*m,3);   % obsolete Sifuentes vector
%v(2*m+1:3*m,1) = -1; v(7*m+1:8*m,2) = -1;   % x,y force balance (g_3,g_4)
%v(1:m,3) = 1; v(5*m+1:6*m,3) = 1;           % mass cons (g_1 + g_2)
%R = randn(2*M,3)/M;                   % randomized version
R = [[ones(M,1);zeros(M,1)],[zeros(M,1);ones(M,1)],[real(p.x);imag(p.x)]]/M;
for i=1:numel(Ns), N = Ns(i);
  s = wormcurve(a,b,Ns(i));
  g = [zeros(2*m,1);jumps(1)+0*L.x;zeros(4*m,1);jumps(2)+0*B.x];  % pres driving
  rhs = [zeros(2*N,1); g]; % driving
  [E,A,Bm,C,Q] = ELSmatrix(s,p,proxyrep,mu,sd,U);
  co = linsolve(E,rhs,lso);
  res(i) = norm(E*co - rhs);
  us(:,i) = evalsol(s,p,proxyrep,mu,sd,U,z,co);   % both cmpts of vel
  Js(:,i) = evalfluxes(s,p,proxyrep,mu,sd,U,co);
  % three non-consistent vectors of inclusion density (Gary)...
  H = [ones(1,N) zeros(1,N);zeros(1,N) ones(1,N);real(s.nx)' imag(s.nx)']';
  %H = randn(2*N,3);                   % randomized version
  Qtilde = Q + (C*H)*R';               % Gary version of Schur
  %if i==1, norm(Q), norm(C*H*R'), svd(Q), svd(Qtilde), end   % sim size norms?
  X = linsolve(Qtilde,C,lso); y = linsolve(Qtilde,g,lso);
  %Bm = Bm + (A*H)*R';     % Gary version of my B corr, preserves nullity 1
  Bm = Bm + (A*H(:,1:2))*R(:,1:2)';  % Alex 
  taut = gmres(@(x) A*x - Bm*(X*x), -Bm*y, [], 1e-14, Ns(i));
  %cond(A - Bm*X)  % 82.7
  xit = y - X*taut;
  taut = taut + H*(R'*xit);         % Gary correction
  cot = [taut;xit];                 % build full soln vector pair
  rest(i) = norm(E*cot - rhs);      % residual back in ELS
  ust(:,i) = evalsol(s,p,proxyrep,mu,sd,U,z,cot);
  Jst(:,i) = evalfluxes(s,p,proxyrep,mu,sd,U,cot);  % Schur flux
  rhsk = knownrhs(src,fsrc,mu,s,U);  % RHS for known 3x3 stokeslet soln...
  cok = linsolve(E,rhsk,lso);   % coeffs for approx to known soln
  %norm(E*cok - rhsk)    % resid for known soln - plot?
  uk = evalsol(s,p,proxyrep,mu,sd,U,z,cok);       % eval this approx
  es(i) = norm(uk-uek);      % err in known case vs exact
end
fprintf('norm X=Qt\\C is %.3g\n',norm(X))
%[U S V] = svd(E); V(:,end)   % show that Nul E = [0;stuff], ie tau unique

disp('u1 N-convergence for ELS, Schur, their diff:')
[us(1,:)',ust(1,:)',us(1,:)'-ust(1,:)']  % 1st cmpt:  u is ELS, ut is Schur
disp('J1 N-convergence for ELS, Schur, their diff:')
[Js(1,:)',Jst(1,:)',Js(1,:)'-Jst(1,:)']
%[Js(2,:)',Jst(2,:)',Js(2,:)'-Jst(2,:)']
figure; %semilogy(Ns,res,'ro-');  % why is resid always around 1e-14?
uerr = sqrt(sum((us-repmat(us(:,end),[1 numel(Ns)])).^2,1)); % 2-norm u errs
semilogy(Ns,uerr,'+-'); hold on;
plot(Ns,abs(Js(1,:)-Js(1,end)),'go-');
uerrt = sqrt(sum((ust-repmat(ust(:,end),[1 numel(Ns)])).^2,1)); % 2-norm ut errs
plot(Ns,uerrt,'+--');
plot(Ns,abs(Jst(1,:)-Jst(1,end)),'go--');
%plot(Ns,res,'r*-'); plot(Ns,rest,'rd-');
plot(Ns,es,'ks-');
legend('u conv ELS','J_1 conv ELS','u conv Schur','J_1 conv Schur','u err vs known');
%legend('u conv ELS','J_1 conv ELS','u conv Schur','J_1 conv Schur','resid ELS','resid Schur','u err vs known');
xlabel('N'); text(70,1e-4,'(c)');
text(140,1e-8, sprintf('$M=%d$,     $m=%d$',M,m),'interpreter','latex');
axis([Ns(1) Ns(end-1) 1e-15 1e-3]);
%set(gcf,'paperposition',[0 0 3.5 3.5]); print -depsc2 figs/stoconvK1.eps
end

if 0, Ms = 10:5:120;    % -------------------- M convergence (not incl Schur)
N = 100; s = wormcurve(a,b,N);  % fixed
g = [zeros(2*m,1);jumps(1)+0*L.x;zeros(4*m,1);jumps(2)+0*B.x];  % pres driving
rhs = [zeros(2*N,1); g]; % driving
Js = nan(2,numel(Ms)); nrms = nan*Ms;
nsings = 6; sings = nan(numel(Ms),nsings);   % save some sing vals
for i=1:numel(Ms)
  p.x = Rp * exp(1i*(1:Ms(i))'/Ms(i)*2*pi); p = setupquad(p);  % reset proxy pts
  E = ELSmatrix(s,p,proxyrep,mu,sd,U);
  S = svd(E);
  sings(i,1) = max(S); sings(i,2:nsings) = S(end-nsings+2:end); % 1st & last few
  co = linsolve(E,rhs,lso);
  nrms(i) = norm(co);
  Js(:,i) = evalfluxes(s,p,proxyrep,mu,sd,U,co);
end
Js(1,:)'
figure; semilogy(Ms,abs(Js(1,:)-Js(1,end)),'b+-'); hold on;
semilogy(Ms,sings(:,2:end),'-','color',.5*[1 1 1]);
semilogy(Ms,nrms,'b.-');
% *** TODO: bring in K=1e3 M-conv data and add to plot as squares
text(15,max(nrms)/10,'(e)');
text(60,max(nrms)/10,sprintf('$N=%d,   m=%d$',N,m),'interpreter','latex');
xlabel('M'); axis([Ms(1) Ms(end-1) 1e-17 max(nrms)]);
%set(gcf,'paperposition',[0 0 3 3]); print -depsc2 figs/stoMconv.eps
end

if 0 % ------------ Schur tests warm-up (see above for convergence)
  N = 140; s = wormshape(a,b,N);    % reset curve
  M = 40; p.x = Rp * exp(1i*(1:M)'/M*2*pi); p = quadr(p);     % reset proxy pts
  rhs = [0*s.x; jumps(1)+0*L.x; 0*L.x; jumps(2)+0*B.x; 0*B.x];   % driving
  [E,~,~,C,Q] = ELSmatrix(s,p,proxyrep,U);
  w = [0*L.w';L.w';0*B.w';B.w'];   % consistency vector in Nul Q^T
  fprintf('norm w^T Q = %.3g\n',norm(w'*Q)) 
  fprintf('norm Q\\C  = %.3g\n',norm(Q\C))
  %svd(Q\C) % just one huge sing val, then gap to a decaying sequence.
  %[U S V] = svd(Q); S = diag(S); r = numel(S); figure; plot(U(:,r),'+-');
end

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

function u = knownsol(U,z,src,fsrc,mu)
% z = targets. U = unit cell struct, src = 1-pt struct, fsrc = source force vec.
% src.nei sets the src grid.  output vel only
[tx ty] = meshgrid(-src.nei:src.nei); trk = tx(:)+1i*ty(:);     % grid
u = srcsum(@StoSLP, trk, [], struct('x',z), src, mu, fsrc);  % dens=fsrc

function rhs = knownrhs(src,fsrc,mu,s,U)
% src is a 1-pt struct with x, w, nx; it will be summed over sec.nei grid,
% unit mag. s is usual target curve (needs normals). fsrc = source force vec.
% The rhs will always be consistent, with f and g nonconstant. Matches knownsol.
[tx ty] = meshgrid(-src.nei:src.nei); trk = tx(:)+1i*ty(:);   % src grid
f = srcsum(@StoSLP,trk,[],s,src,mu, fsrc);  % direct vel data sum to curve
U.nei = src.nei;                  % sets up C correct for src grid
g = Cblock(src,U,mu,[1 0]) * fsrc;               % sd sets SLP only
rhs = [f;g];

function [zz ii] = extgrid(gx,gy,s,U)  % grid points and indices outside Omegas
% given gx,gy 1d x,y grids, s=curve segment, U = unit cell struct, return
% zz = col list of grid pts as C-#s, and ii = logical index array if outside
% Omega and all images
[xx yy] = meshgrid(gx,gy); zz = (xx(:)+1i*yy(:)); clear xx yy;
ii = true(size(zz));               % indices outside Omega or its copies
for i=-1:1, for j=-1:1, si = s.inside(zz+U.e1*i+U.e2*j); ii(si)=false; end, end
