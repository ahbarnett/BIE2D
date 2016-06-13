% Basic solver for Stokes interior velocity (Dirichlet) BVP on smooth curve.
% Also serves to test the native SLP and DLP matrix and eval routines.
% Barnett 6/13/16, munging LapintDirBVP.m with intvelBVP.m from Feb 2014.
%
% References: see Sec 2.3 of this:
%  [HW]    "Boundary Integral Equations", G. C. Hsiao and W. L. Wendland
%          (Springer, 2008).

clear;
mu = 0.7;              % viscosity
a = .3; w = 5;         % smooth wobbly radial shape params
N = 200; s = wobblycurve(a,w,N);

p = []; p.x = [.2+.3i; .1-.4i]; p.nx = exp([1.9i;-0.6i]);  % 2 test pts w/ dirs
%figure; showsegment({s p})

% known Stokes solution: Poisseuil pipe flow
ue = @(x) [imag(x).^2;0*x]; pe = @(x) 2*mu*real(x);   % vel, pres
Te = @(x,n) 2*mu*[-real(x).*real(n)+imag(x).*imag(n); -real(x).*imag(n)+imag(x).*real(n)];   % trac (n=targ nor)

% double-layer representation...
A = -eye(2*N)/2 + StoDLPmat(s,s);       % Nystrom discretized integral operator
%S = svd(A); fprintf('last few singular values of D-1/2:\n'); S(end-4:end)
A(:,1) = A(:,1) + [real(s.nx);imag(s.nx)];   % rank-1 perturbation kills nul A
%S = svd(A); fprintf('last few singular values of A:\n'); S(end-4:end)
rhs = ue(s.x);                          % velocity data
tau = A \ rhs;
fprintf('resid norm %.3g,  density norm %.3g\n',norm(rhs-A*tau),norm(tau))

[up pp Tp] = StoDLPeval(p,s,tau,mu);    % vel, pres, trac, at test pts
fprintf('D rep:\tnative u error @ test pts:\t\t%.3g\n',max(abs(up-ue(p.x))))
fprintf('\tnative p diff error btw test pts: \t%.3g\n',diff(pp)-diff(pe(p.x)))
poff = mean(pp - pe(p.x));                   % get observed pres offset in rep
Tp = Tp - poff*[-real(p.nx);-imag(p.nx)];      % fix traction's pres offset
fprintf('\tnative traction err @ test pts: \t%.3g\n',max(abs(Tp-Te(p.x,p.nx))))

return

% *** finish the plots & S test....

fnp = fx(p.x)*real(p.nx) + fy(p.x)*imag(p.nx);   % exact targ n-deriv 


% plot solution on grid & soln errors...
nx = 200; gx = max(abs(s.x))*linspace(-1,1,nx);
[xx yy] = meshgrid(gx); g.x = xx(:)+1i*yy(:);
ug = LapDLPeval(g,s,tau);       % u on grid, expensive bit for now
fg = f(g.x);                    % known soln on grid
ug = reshape(ug,[nx nx]); fg = reshape(fg,[nx nx]);  % shape arrays for plot
figure; set(gcf,'name', 'DLP native evaluation');
tsubplot(1,2,1); imagesc(gx,gx,ug); showsegment(s);
caxis([-1 2]); colorbar; axis tight; title('u');
tsubplot(1,2,2); imagesc(gx,gx,log10(abs(ug-fg))); showsegment(s);
caxis([-16 0]); colorbar; axis tight; title('log_{10} error u');
% basic BVP done

% instead use close-evaluation scheme...
ii = s.inside(g.x); g.x = g.x(ii); ug = nan*ug;  % eval only at interior pts
ug(ii) = LapDLPeval_closeglobal(g,s,tau,'i');
figure; set(gcf,'name', 'DLP close evaluation');
tsubplot(1,2,1); imagesc(gx,gx,ug); showsegment(s);
caxis([-1 2]); colorbar; axis tight; title('u');
tsubplot(1,2,2); imagesc(gx,gx,log10(abs(ug-fg))); showsegment(s);
caxis('auto'); colorbar; axis tight; title('log_{10} error u');

% mixed double plus single rep... (not helpful---cond(A) worse---but tests SLP)
A = A + LapSLPmat(s,s);
tau = A \ rhs;
fprintf('resid norm %.3g,  density norm %.3g\n',norm(rhs-A*tau),norm(tau))
[up unp] = LapDLPeval(p,s,tau); [vp vnp] = LapSLPeval(p,s,tau);
up = up+vp; unp = unp+vnp;
fprintf('D+S rep: native u and un errors @ test pt: \t%.3g\t%.3g \n',up-f(p.x),unp-fnp)

