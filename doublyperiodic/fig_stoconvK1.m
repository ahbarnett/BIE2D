function fig_stoconvK1
% self-contained code for Stokes no-slip periodic, convergence and soln plots
% Single inclusion (K=1), no close eval.  Adapted from fig_lapconvK1.m,
% except SLP-type proxies hard-wired, and close eval done.  Barnett 6/7/16

warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square
mu = 0.7;   % overall viscosity const for Stokes PDE
sd = [1 1];  % inclusion representation prefactors for SLP+DLP
jumps = [1 0]; %[0 1]; % pressure jumps across R-L and T-B (not for known soln)

U.e1 = 1; U.e2 = 1i; U.nei = 1;  % unit cell; nei=1 for 3x3 scheme
m = 22; [U L R B T] = walls(U,m);
Rp = 1.4; M = 80;        % proxy params
p.x = Rp * exp(1i*(1:M)'/M*2*pi); p = quadr(p); % proxy pts

a = 0.7; b = 0.15;   % worm params, spills horizontally out of any unit cell
uex = [0.016778793238;0.005152952237]; flux1ex = nan;  % a=.7,b=.15;  to 1e-12
% known soln stokeslet location: must be deep inside Omega (careful)
src.x = 0.42+.23i; src.w = 1; src.nx = 1+0i;

% -------------------------- single soln and plot
N =150; s = wormshape(a,b,N);
% obstacle no-slip & pressure-drop driving...
rhs = [zeros(2*N,1); zeros(2*m,1);jumps(1)+0*L.x;zeros(4*m,1);jumps(2)+0*B.x];
tic
E = ELSmatrix(s,p,@SLPmatrix,U,mu,sd);              % fill
co = linsolve(E,rhs,lso);                           % direct bkw stable solve
toc
%S = svd(E); disp('last few sing vals of E:'), S(end-5:end) % dim Nul E = 1
fprintf('resid norm = %.3g\n',norm(E*co - rhs))
sig = co(1:2*N); psi = co(N+1:end);
fprintf('density norm = %.3g, proxy norm = %.3g\n',norm(sig), norm(psi))
fprintf('body force + jumps = (%.3g,%.3g)  should vanish\n',s.w'*sig(1:N)+abs(U.e1)*jumps(1),s.w'*sig(N+1:end)+abs(U.e2)*jumps(2))
z = .1+.4i;
u = evalvel(s,p,mu,sd,U,z,co);
fprintf('u at pt = (%.16g,%.16g)  \t(est abs err: %.3g)\n',u(1),u(2),norm(u-uex))
%tic; J = evalfluxes(s,p,@SLPmatrix,U,co); toc
%fprintf('fluxes (%.16g,%.16g)   (est abs err: %.3g)\n',J(1),J(2),J(1)-flux1ex)
if 1   % figure
  nx = 101; gx = 0.5*((0:nx-1)/(nx-1)*2-1); ng = nx^2;  % fine pres eval grid
  gy = gx; [xx yy] = meshgrid(gx,gy); zz = (xx(:)+1i*yy(:)); clear xx yy;
  tic; pg = evalpres(s,p,mu,sd,U,zz,co,0); toc
  pg = reshape(pg,[nx nx]) - pg(ceil(ng/2));   % offset by pres @ 0
  for i=-1:1, for j=-1:1, pg(s.inside(zz+U.e1*i+U.e2*j)) = nan; end, end % tidy
  figure; contourf(gx,gy,pg,[-.6:.1:.6]); colormap(jet(256)); hold on;  % pres
  nx = 26; gx = 0.5*((0:nx-1)/(nx-1)*2-1); ng = nx^2;  % coarse vel eval grid
  gy = gx; [xx yy] = meshgrid(gx,gy); zz = (xx(:)+1i*yy(:)); clear xx yy;
  in = 0*zz;                          % indices inside (not to evaluate)
  for i=-1:1, for j=-1:1, si = s.inside(zz+U.e1*i+U.e2*j); in(si) = 1; end, end
  ug = nan(2*ng,1);
  tic; ug([~in;~in]) = evalvel(s,p,mu,sd,U,zz(~in),co,1); toc
  u1 = reshape(ug(1:ng),[nx nx]); u2 = reshape(ug(ng+1:end),[nx nx]); % u cmpts
  quiver(gx,gy,u1,u2,2.0,'k-');       % show vector field
  n=1; for i=-n:n, for j=-n:n, plot([s.x;s.x(1)]+j+1i*i,'-'); end,end % curves
  plot([L.x;R.x;T.x;B.x],'b.');  % wall nodes
  axis xy equal off; plot(z,'k.'); axis([0 1 0 1]-0.5);
  %for i=-1:1, for j=-1:1, plot(src.x+U.e1*i+U.e2*j,'k*'); end, end % known
  text(0,0,'$\Omega$','interpreter','latex');
  %text(-1.2,1.3,'(a)'); %,'fontsize',14);
  %set(gcf,'paperposition',[0 0 4 4]); print -depsc2 stosolK1.eps
end

keyboard; return







if 1, Ns = 30:10:230;   % ------------------------  N-convergence
us = nan*Ns; res = us; rest = us; ust = us; es = us;
Js = nan(2,numel(Ns)); Jst = Js;
uek = knownsol(U,z,src);   % known dipole 3x3 grid soln, fixed, ignores jumps
v = [0*L.w';1+0*L.w';0*B.w';1+0*B.w']; d = ones(M,1)/M; % Sifuentes vectors
for i=1:numel(Ns)
  s = wormshape(a,b,Ns(i));
  g = [jumps(1)+0*L.x; 0*L.x; jumps(2)+0*B.x; 0*B.x]; rhs = [0*s.x; g]; %driving
  [E,A,Bm,C,Q] = ELSmatrix(s,p,proxyrep,U);
  co = linsolve(E,rhs,lso);
  res(i) = norm(E*co - rhs);
  u = evalsol(s,p,proxyrep,U,z,co);
  us(i) = u(2)-u(1);
  Js(:,i) = evalfluxes(s,p,proxyrep,U,co);
  Qtilde = Q + v*d';  % Schur stuff...
  %if i==1, norm(Q), norm(v*d'), svd(Q), svd(Qtilde), end   % sim size norms?
  QtdagC = linsolve(Qtilde,C,lso); Qtdagg = linsolve(Qtilde,g,lso);
  %Qtdag = pinv(Qtilde); QtdagC = Qtdag*C; Qtdagg = Qtdag*g;  % loses 7 digits
  taut = gmres(@(x) A*x - Bm*(QtdagC*x), -Bm*Qtdagg, [], 1e-14, Ns(i));
  %taut = linsolve(A - Bm*QtdagC,-Bm*Qtdagg,lso);  % direct soln
  %cond(A - Bm*QtdagC)  % 8.4
  xit = Qtdagg - QtdagC*taut; cot = [taut;xit];  % build full soln vector
  rest(i) = norm(E*cot - rhs);      % residual back in ELS
  u = evalsol(s,p,proxyrep,U,z,cot);
  ust(i) = u(2)-u(1);
  Jst(:,i) = evalfluxes(s,p,proxyrep,U,cot);  % Schur flux
  rhsk = knownrhs(src,s,U);     % set up RHS for known 3x3 unit source soln...
  cok = linsolve(E,rhsk,lso);   % coeffs for approx to known soln
  %norm(E*cok - rhsk)    % resid for known soln - plot?
  uk = evalsol(s,p,proxyrep,U,z,cok);       % eval this approx
  es(i) = uk(2)-uk(1)-(uek(2)-uek(1));      % err vs known diff btw test pts
end
fprintf('norm Qt\\C  = %.3g\n',norm(Qtilde\C))
[us',ust',us'-ust']
[Js(1,:)',Jst(1,:)',Js(1,:)'-Jst(1,:)']
[Js(2,:)',Jst(2,:)',Js(2,:)'-Jst(2,:)']
figure; %semilogy(Ns,res,'ro-');  % why is resid always around 1e-14?
semilogy(Ns,abs(us-us(i)),'+-'); hold on;
plot(Ns,abs(Js(1,:)-Js(1,end)),'go-');
plot(Ns,abs(ust-ust(i)),'+--');
plot(Ns,abs(Jst(1,:)-Jst(1,end)),'go--');
%plot(Ns,rest,'rd-');
plot(Ns,abs(es),'ks-');
%legend('u convergence','J_1 convergence','u err vs known');
legend('u conv ELS','J_1 conv ELS','u conv Schur','J_1 conv Schur','u err vs known');
xlabel('N'); text(40,1e-4,'(c)');
text(140,1e-8, sprintf('$M=%d$,     $m=%d$',M,m),'interpreter','latex');
axis([Ns(1) Ns(end-1) 1e-15 1e-3]);
%set(gcf,'paperposition',[0 0 3.5 3.5]); print -depsc2 lapconvK1.eps
end

%keyboard

if 0, Ms = 10:5:120;    % -------------------- M convergence, incl Schur
N = 100; s = wormshape(a,b,N);  % fixed
rhs = [0*s.x; jumps(1)+0*L.x; 0*L.x; jumps(2)+0*B.x; 0*B.x];   % driving
Js = nan(2,numel(Ms)); nrms = nan*Ms;
nsings = 6; sings = nan(numel(Ms),nsings);   % save some sing vals
for i=1:numel(Ms)
  p.x = Rp * exp(1i*(1:Ms(i))'/Ms(i)*2*pi); p = quadr(p);     % reset proxy pts
  E = ELSmatrix(s,p,proxyrep,U);
  S = svd(E);
  sings(i,1) = max(S); sings(i,2:nsings) = S(end-nsings+2:end); % 1st & last few
  co = linsolve(E,rhs,lso);
  nrms(i) = norm(co);
  Js(:,i) = evalfluxes(s,p,proxyrep,U,co);
end
Js(1,:)'
figure; semilogy(Ms,abs(Js(1,:)-Js(1,end)),'b+-'); hold on;
semilogy(Ms,sings(:,2:end),'-','color',.5*[1 1 1]);
semilogy(Ms,nrms,'b.-');
% *** TODO: bring in K=1e3 M-conv data and add to plot as squares
text(15,max(nrms)/10,'(e)');
text(60,max(nrms)/10,sprintf('$N=%d,   m=%d$',N,m),'interpreter','latex');
xlabel('M'); axis([Ms(1) Ms(end-1) 1e-17 max(nrms)]);
set(gcf,'paperposition',[0 0 3 3]); print -depsc2 lapMconv.eps
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


function [U L R B T] = walls(U,M)  % setup walls with M nodes for unit square U
[x w] = gauss(M); w=w/2;
L.x = (-U.e1 + U.e2*x)/2; L.nx = (-1i*U.e2)/abs(U.e2) + 0*L.x; L.w=w*abs(U.e2);
R = L; R.x = L.x + U.e1;
B.x = (-U.e2 + U.e1*x)/2; B.nx = (1i*U.e1)/abs(U.e1) + 0*B.x; B.w=w*abs(U.e1);
T = B; T.x = B.x + U.e2;
U.L = L; U.T = T; U.R = R; U.B = B;

function u = evalvel(s,p,mu,sd,U,z,co,close)     % evaluate velocity field
% z = list of targets as C values. p = proxy struct, s = source curve struct
% co = full coeff vec, U = unit cell struct, mu = viscosity
% sd = prefactors for source rep, close enables special close quadr. 6/7/16
if nargin<8, close=0; end
N = numel(s.x); nei = U.nei;
sig = co(1:2*N); psi = co(2*N+1:end); % split up solution coeffs (psi=proxy)
u = SLPmatrix(struct('x',z),p,mu) * psi;      % init u w/ proxies (always far)
if ~close                             % naive quadr rule
  z = struct('x',z);
  for i=-nei:nei, for j=-nei:nei, a = U.e1*i+U.e2*j;    % src trans
      u = u + sd(1)*(SLPmatrix(z,s,mu,a) * sig) + sd(2)*(DLPmatrix(z,s,mu,a) * sig);  % naive matvecs, vel due to completed curve rep
    end,end
else                                  % use exterior close eval for all pts
  sigc = sig(1:N) + 1i*sig(N+1:end);  % pack as C-format density
  for i=-nei:nei, for j=-nei:nei, a = U.e1*i+U.e2*j; % src trans (done on targ)
      uij = (sd(1)/mu)*StokesScloseeval(z-a,s,sigc,'e') + sd(2)*StokesDcloseeval(z-a,s,sigc,'e');   % note mu factor only for SLP defn
      u = u + [real(uij);imag(uij)];  % unpack C-format vel
    end,end
end
 
function pr = evalpres(s,p,mu,sd,U,z,co,close)
% same as evalvel but evaluates pressure. 6/7/16. *** TODO CLOSE EVAL
if nargin<8, close=0; end
z = struct('x',z); N = numel(s.x); nei = U.nei;
sig = co(1:2*N); psi = co(2*N+1:end);    % split up solution coeffs (psi=proxy)
pr = SLPpresmatrix(z,p,mu) * psi;        % init w/ proxies
for i=-nei:nei, for j=-nei:nei, a = U.e1*i+U.e2*j;    % src trans
    pr = pr + sd(1)*(SLPpresmatrix(z,s,mu,a) * sig) + sd(2)*(DLPpresmatrix(z,s,mu,a) * sig);  % pressure due to completed curve rep
  end,end

function J = evalfluxes(s,p,proxyrep,U,co)   % ***** MAKE STOKES
% inputs as in evalsol. Uses Gary-inspired bdry of 3x3 block far-field method
if U.nei~=1, warning('U.neu must equal 1'); end
w = U.L.w; if norm(w-U.B.w)>1e-14, error('L and B must have same weights'); end
m = numel(w);
N = numel(s.x); sig = co(1:N); psi = co(N+1:end);
t.x = [U.L.x;U.B.x]; t.nx = [U.L.nx;U.B.nx];  % 2-wall target
[~,Tn] = proxyrep(t, p); u = Tn * psi;   % proxy contrib to un only
J = sum(reshape(u,[m 2])'.*([1;1]*w),2);         % .... do its quadr on L,B 
for i=-1:1      % set up big loop of 12 walls, just their nodes
  x{i+2} = U.L.x - U.e1 +i*U.e2; x{i+5} = x{i+2} + 3*U.e1;
  x{i+8} = U.B.x +i*U.e1 -U.e2; x{i+11} = x{i+8} + 3*U.e2;
end
t.x = vertcat(x{:}); t.nx = [repmat(U.L.nx,[6 1]); repmat(U.B.nx,[6 1])];
[~,Tn] = SLPmatrix(t, s); u = Tn * sig;   % central density only
amts = [0 0 0 3 3 3 -1 -2 -3 1 2 3; -1 -2 -3 1 2 3 0 0 0 3 3 3];  % wall wgts
J = J + sum(repmat(u',[2 1]).*kron(amts,w),2);   % weight each wall

function s = wormshape(a,b,N)    % define curve and quadr.   a,b shape params
s.t = (1:N)'/N*2*pi; s.x = a*cos(s.t)+b*1i*sin(s.t);
s.x = s.x + 0.3i*sin(2*real(s.x));   % spec diff can limit acc in this shape
s.inside = @(z) (real(z)/a).^2+((imag(z)-0.3*sin(2*real(z)))/b).^2 < 1; % yuk
s = quadr(s);  % uses spectral differentiation, not as acc as analytic kappa(s)

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

function h = showseg(s, U)    % plot segments & possibly its 2D nei images
if nargin<2, U.nei = 0; U.e1 = 0; U.e2 = 0; end % dummy uc
if iscell(s), for i=1:numel(s), showseg(s{i},U); end, return, end % multi segs
for i=-U.nei:U.nei, for j=-U.nei:U.nei
  plot(s.x+U.e1*i+U.e2*j, 'b.-'); axis equal xy;
  l=0.05; hold on; plot([s.x+U.e1*i+U.e2*j, s.x+l*s.nx+U.e1*i+U.e2*j].', 'k-');
end, end

function [E A B C Q] = ELSmatrix(s,p,proxyrep,U,mu,sd)
% builds matrix blocks for Stokes extended linear system, S+D rep w/ Kress
nei = U.nei; N = numel(s.x);
A = (sd(2)/2)*eye(2*N);    % exterior S+D self vel matrix, starting w/ jump rel
for i=-nei:nei, for j=-nei:nei  % direct sum
    A = A + sd(1)*SLPmatrix(s,s,mu,U.e1*i+U.e2*j) + sd(2)*DLPmatrix(s,s,mu,U.e1*i+U.e2*j);
  end,end
B = proxyrep(s,p,mu);     % map from proxy density to vel on curve
C = Cblock(s,U,mu,sd);
[QL QLt] = proxyrep(U.L,p,mu); [QR QRt] = proxyrep(U.R,p,mu); % vel, traction
[QB QBt] = proxyrep(U.B,p,mu); [QT QTt] = proxyrep(U.T,p,mu);
Q = [QR-QL; QRt-QLt; QT-QB; QTt-QBt];
E = [A B; C Q];

function C = Cblock(s,U,mu,sd)     % fill C from source curve s to U walls
% sd controls prefactors on SLP & DLP for Stokes rep on the obstacle curve
nei = U.nei; N = numel(s.x); m = numel(U.L.x);
CL = zeros(2*m,2*N); CLn = CL;   % by "n" we of course mean traction
for i=-nei:nei
  [Ci Cni] = SLPmatrix(U.L,s,mu,nei*U.e1+i*U.e2);
  CL = CL+sd(1)*Ci; CLn=CLn+sd(1)*Cni;
  [Ci Cni] = DLPmatrix(U.L,s,mu,nei*U.e1+i*U.e2);
  CL = CL+sd(2)*Ci; CLn=CLn+sd(2)*Cni;
end
CR = 0*CL; CRn = CR;
for i=-nei:nei
  [Ci Cni] = SLPmatrix(U.R,s,mu,-nei*U.e1+i*U.e2);
  CR = CR+sd(1)*Ci; CRn=CRn+sd(1)*Cni;
  [Ci Cni] = DLPmatrix(U.R,s,mu,-nei*U.e1+i*U.e2);
  CR = CR+sd(2)*Ci; CRn=CRn+sd(2)*Cni;
end
CB = 0*CL; CBn = CB;
for i=-nei:nei
  [Ci Cni] = SLPmatrix(U.B,s,mu,nei*U.e2+i*U.e1);
  CB = CB+sd(1)*Ci; CBn=CBn+sd(1)*Cni;
  [Ci Cni] = DLPmatrix(U.B,s,mu,nei*U.e2+i*U.e1);
  CB = CB+sd(2)*Ci; CBn=CBn+sd(2)*Cni;
end
CT = 0*CL; CTn = CT;
for i=-nei:nei
  [Ci Cni] = SLPmatrix(U.T,s,mu,-nei*U.e2+i*U.e1);
  CT = CT+sd(1)*Ci; CTn=CTn+sd(1)*Cni;
  [Ci Cni] = DLPmatrix(U.T,s,mu,-nei*U.e2+i*U.e1);
  CT = CT+sd(2)*Ci; CTn=CTn+sd(2)*Cni;
end
C = [CR-CL; CRn-CLn; CT-CB; CTn-CBn];

function u = knownsol(U,z,src)
% z = targets. U = unit cell struct, src = charge loc.  output potential only
% Note use of dipole (monopole would work, but dipole closer to application)
nei = U.nei;
u = 0*z;
for i=-nei:nei, for j=-nei:nei
    u = u + DLPmatrix(struct('x',z), src, U.e1*i+U.e2*j);  % pot due to DLP
  end,end

function rhs = knownrhs(src,s,U)
% src is a dipole source with x, w, nx; it will be summed over 3x3, unit mag.
% s is usual target curve (needs normals).
% The rhs will always be consistent, with f and g nonconstant. Matches knownsol.
nei = U.nei;
A = zeros(numel(s.x),numel(src.x));
for i=-nei:nei, for j=-nei:nei        % direct sum
    [~,Aij] = DLPmatrix(s,src,U.e1*i+U.e2*j); A = A + Aij;
  end,end
f = A;
g = Cblock(src,U,@DLPmatrix);
rhs = [f;g];

function [A,T] = SLPmatrix(t,s,mu,a) % single-layer 2D Stokes kernel vel matrix
% Returns 2N-by-2N matrix from src force vector to 2 flow component
% t = target seg (x cols), s = src seg, a = optional translation of src seg.
% Kress log-singularity quadrature for self-interaction.
% 2nd output is traction matrix, needs t.nx normal (C-#); no self-int yet.
if nargin<4, a = 0; end  % complex number interpreted as x+iy
N = numel(s.x); M = numel(t.x);
r = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
irr = 1./(conj(r).*r);    % 1/r^2
d1 = real(r); d2 = imag(r);
logir = -log(abs(r));  % log(1/r) diag block
c = 1/(4*pi*mu);       % factor from Hsiao-Wendland book, Ladyzhenskaya
if numel(s.x)==numel(t.x) && max(abs(s.x+a-t.x))<1e-14   % self via Kress
  S = logir + circulant(0.5*log(4*sin(pi*(0:N-1)/N).^2)); % peri log
  S(diagind(S)) = -log(s.sp);                       % diagonal limit
  m = 1:N/2-1; Rjn = ifft([0 1./m 2/N 1./m(end:-1:1)])/2;  % Kress Rj(N/2)/4pi
  S = S/N + circulant(Rjn); % includes SLP prefac 1/2pi. Kress peri log matrix L
  S = S .* repmat(s.sp.',[N 1]);  % include speed factors (not 2pi/N weights)
  A = (1/2/mu) * kron(eye(2),S);       % prefactor & diagonal blocks
  t1 = real(s.tang); t2 = imag(s.tang);  % now do r tensor r part...
  A11 =  d1.^2.*irr; A11(diagind(A11)) = t1.^2;     % diagonal limits
  A12 =  d1.*d2.*irr; A12(diagind(A12)) = t1.*t2;
  A22 =  d2.^2.*irr; A22(diagind(A22)) = t2.^2;
  A = A + c*[A11 A12; A12 A22].*repmat([s.w(:)' s.w(:)'], [2*M 1]); % pref & wei
else                     % distinct src and targ
  A12 = d1.*d2.*irr;     % off diag vel block
  A = c*[logir + d1.^2.*irr, A12;                         % u_x
         A12,                logir + d2.^2.*irr];         % u_y
  A = A .* repmat([s.w(:)' s.w(:)'], [2*M 1]);            % quadr wei
end
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
A = (1/2/pi) * [d1.*irr, d2.*irr];                   % pressure
A = A .* repmat([s.w(:)' s.w(:)'], [M 1]);           % quadr wei

function [A,T] = DLPmatrix(t,s,mu,a) % double-layer 2D Stokes vel kernel matrix
% Returns 2N-by-2N matrix from src force vector to 2 flow components.
% If detects self-int, does correct diagonal limit, no jump condition (PV int).
% t = target seg (x cols), s = src seg, a = optional transl of src seg.
% 2nd output is optional traction matrix, without self-eval. Barnett 3/2/14
if nargin<4, a = 0; end  % complex number interpreted as x+iy
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
