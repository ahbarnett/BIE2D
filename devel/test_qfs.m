% Quadrature by Fundamental Solutions (QFS) for
% Laplace interior Dirichlet BVP on smooth curve.
% Barnett 2/12/19. Hacked from solvers/LapintDirBVP.m

run('../bie2dsetup')
clear; a = .3; w = 5;         % smooth wobbly radial shape params
N = 480; s = wobblycurve(1,a,w,N);

p = []; p.x = 1+1.1i; p.nx = exp(1.9i);   % ext test pt including a direction
% for this shape & rhs, N=240 gets 1e-8 fou decay @ N/2, 480 gets full 1e-15.
close all; %figure; showsegment({s p})

%z0 = -.1+.1i; f0 = exp(1i*4.0);   % int src pt
z0 = .2+0.7i; f0 = exp(1i*4.0);   % int src pt
fholom = @(z) f0./(z-z0);    % holomorphic ext soln
fpholom = @(z) -f0./(z-z0).^2;   % its complex deriv
f = @(z) real(fholom(z));           % use real part as known Laplace soln
fx = @(z) real(fpholom(z)); fy = @(z) -imag(fpholom(z)); % partials, note sign!

% ------- u = double-layer representation...  (ext BVP has rank-1 deficient)
L = ones(N,1)*s.w'/sum(s.w);       % rank-1 pert for well-cond (added Aug 2019)
eta = 0.0;
K = eye(N)/2 + LapDLP(s,s);        % on-surf exterior self-eval mat
if eta~=0, K = K + LapSLP(s,s); end
A = L + K;    % Nystrom discretized integral operator
rhs = f(s.x);
tau = A \ rhs;
fprintf('Lap ext Dir BVP: cond = %.3g,  resid norm %.3g,  density norm %.3g\n',cond(A),norm(rhs-A*tau),norm(tau))
[up unp] = LapDLP(p,s,tau);        % potential and targ n-deriv at test pt
fnp = fx(p.x)*real(p.nx) + fy(p.x)*imag(p.nx);   % exact targ n-deriv
fprintf('D rep: native u and un errors @ test pt: \t%.3g\t%.3g \n',up-f(p.x),unp-fnp)

% plot solution on grid & soln errors...
x0 = 1.2; nx = 200; gx = max(abs(s.x))*linspace(-x0,x0,nx);
[xx yy] = meshgrid(gx); zz = xx+1i*yy; ii = ~s.inside(zz);  % outside
g.x = zz(ii(:));
ug = nan*xx; tic; ug(ii) = LapDLP(g,s,tau); toc     % u on grid, expensive
fg = nan*xx; fg(ii) = f(g.x);                       % known soln on grid
figure(1); set(gcf,'name', 'Lap DLP native eval');
tsubplot(1,3,1); imagesc(gx,gx,fg); showsegment(s);
colorbar; axis tight; title('u exact');
hold on; plot(real(z0),imag(z0),'w*');
tsubplot(1,3,2); imagesc(gx,gx,log10(abs(ug-fg))); showsegment(s);
caxis([-16 0]); colorbar; axis tight; title('Native rule: log_{10} error u');
subplot(1,3,3); semilogy(abs(ifft(tau))); axis tight; ylabel('abs fou coeffs of \tau');
% basic BVP done ------

% precompute QFS for the rigid obstacle...
upsamp = 3.0;    % QFS param, eg 3.0  - make diff for on/off surf.
Nf = ceil(upsamp*N); tf = (0:Nf-1)'*(2*pi/Nf);
onsurf = 0;   % 0: go off-surf, 1: use on-surf rule, in constructing QFS
if onsurf
  e = wobblycurve(1,a,w,Nf);   % same curve, maybe upsampled
  C = LapDLP(e,e) + eye(Nf)/2;   % on-surf eval matrix, w/ ext jump
  Ne = Nf; sf = e;             % # on targ curve, here same
else                   % use nearby-curve
  Ne = N;               % don't upsamp # ext targs
  d0 = 5.0;      % how close in upsamp h-units
  e = []; e.x = s.Z(s.t - 2i*pi/Nf*d0);  % ext targ curve
  sf = wobblycurve(1,a,w,Nf);   % fine src curve
  tauf = perispecinterp(tau,Nf); ue = LapDLP(e,sf,tauf);  % test it
  fprintf('worst upsamp native err on e curve: %.3g\n',max(abs(ue-f(e.x))))
  C = LapDLP(e,sf);
end
I = nan(Nf,N); for i=1:N, v = zeros(N,1); v(i) = 1;  % build interp matrix
  I(:,i) = perispecinterp(v,Nf); end
C = C*I;    % now C maps N-pt dens to ext  curve

% set up FS..
%Nr = ceil(1.3*N);       % needed >N to get good QFS acc matching Cauchy.
Nr = ceil(1.0*N);
tr = (0:Nr-1)'*(2*pi/Nr);
r = []; dr = -5.0; r.x = s.Z(tr - 2i*pi/Nr*dr);  % int FS curve (hr-dist = dr)
r.w = 2*pi/Nr;  % dummy weights for pt src FS
figure(2); showsegment({sf,e,r}); hold on; plot(real(z0),imag(z0),'r*');
B = LapSLP(e,r);

factor = 0;   % 0: single pinv matrix, 1: two factors (needed for >7 digits)
if ~factor
  E = B\C;   % eval matrix from density on sf to FS source strengths, too big
  %[L,U] = lu(B); E = U\(L\C);
else
  eps = 1e-15; [U,S,V] = svd(B);
  ra = sum(diag(S)>eps); S = diag(S); S = S(1:ra); iS = 1./S; % rank
  Bi2 = V(:,1:ra)*diag(iS); Bi1 = U(:,1:ra); Bi1 = Bi1';    % two factors of B^+
  % check soln of lin sys Bx = y...
  %y = f(e.x); x1 = B\y; [norm(x1), norm(B*x1-y)]
  %x2 = Bi2*(Bi1*y); [norm(x2), norm(B*x2-y)]
  Ei1 = Bi1*C;        % now the pair (Ei1, Bi2) are all we need
end
  
% Test close-evaluation schemes...
tic; ug(ii) = LapDLP_closeglobal(g,s,tau,'e'); toc
fprintf('max grid DLP pot close eval soln err:\t%.3g\n',max(abs(ug(:)-fg(:))))
figure(3); set(gcf,'name', 'Lap DLP close eval');
tsubplot(1,2,1); imagesc(gx,gx,log10(abs(ug-fg))); showsegment(s);
caxis('auto'); colorbar; axis tight; title('Cauchy-comp: log_{10} err u');
fprintf('max Cauchy-comp err on grid: %.3g\n',max(abs(ug(:)-fg(:))))

if ~factor
  sig = E*tau;          % apply QFS as single pinv - loses digits
  fprintf('||E||_2 = %.3g\n',norm(E))
else
  sig = Bi2*(Ei1*tau);  % apply QFS: two factor
  fprintf('||Ei1||_2 = %.3g,   ||Bi2||_2 = %.3g\n',norm(Ei1),norm(Bi2))
end
fprintf('QFS: norm(sig) = %.3g\n',norm(sig))
tic; ug(ii) = LapSLP(g,r,sig); toc   % eval the FS
tsubplot(1,2,2); imagesc(gx,gx,log10(abs(ug-fg))); showsegment(s);
caxis('auto'); colorbar; axis tight; title('new QFS: log_{10} err u');
fprintf('max QFS err on grid: %.3g\n',max(abs(ug(:)-fg(:))))
