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
N = 300; s = wobblycurve(a,w,N);

p = []; p.x = [.2+.3i; .1-.4i]; p.nx = exp([1.9i;-0.6i]);  % 2 test pts w/ dirs
%figure; showsegment({s p})

% known Stokes solution: Poisseuil pipe flow
ue = @(x) [imag(x).^2;0*x]; pe = @(x) 2*mu*real(x);   % vel, pres
Te = @(x,n) 2*mu*[-real(x).*real(n)+imag(x).*imag(n); -real(x).*imag(n)+imag(x).*real(n)];   % trac (n=targ nor)

% double-layer representation...
A = -eye(2*N)/2 + StoDLP(s,s,mu);       % Nystrom discretized integral operator
%S = svd(A); fprintf('last few singular values of D-1/2:\n'); S(end-4:end)
A(:,1) = A(:,1) + [real(s.nx);imag(s.nx)];   % rank-1 perturbation kills nul A
%S = svd(A); fprintf('last few singular values of A:\n'); S(end-4:end)
rhs = ue(s.x);                          % velocity data
tau = A \ rhs;
fprintf('Sto int Dir BVP resid norm %.3g,  density norm %.3g\n',norm(rhs-A*tau),norm(tau))

[up pp Tp] = StoDLP(p,s,mu,tau);        % vel, pres, trac, at test pts
fprintf('D rep:\tnative u error @ test pts:\t\t%.3g\n',max(abs(up-ue(p.x))))
fprintf('\tnative p diff error btw test pts: \t%.3g\n',diff(pp)-diff(pe(p.x)))
poff = mean(pp - pe(p.x));                   % get observed pres offset in rep
Tp = Tp - poff*[-real(p.nx);-imag(p.nx)];    % fix traction's pres offset
fprintf('\tnative traction err @ test pts: \t%.3g\n',max(abs(Tp-Te(p.x,p.nx))))

% plot solution on grid & soln errors...
nx = 100; gx = max(abs(s.x))*linspace(-1,1,nx);
[xx yy] = meshgrid(gx); g.x = xx(:)+1i*yy(:);
tic, [ug pg] = StoDLP(g,s,mu,tau); toc    % u, p on grid, expensive bit for now
ueg = ue(g.x); peg = pe(g.x);             % known soln on grid
ug = reshape(ug,[nx nx 2]); ueg = reshape(ueg,[nx nx 2]);
pg = reshape(pg,[nx nx]); peg = reshape(peg,[nx nx]) + poff;   % fix p offset
figure; set(gcf,'name', 'Sto DLP native u,p eval');
spg = abs(ug(:,:,1)+1i*ug(:,:,2));   % flow speed on grid
tsubplot(2,2,1); imagesc(gx,gx,spg); showsegment(s);
caxis([0 1.5]); colorbar; axis tight; title('|u|');
uerrg = abs(ug(:,:,1)+1i*ug(:,:,2)-ueg(:,:,1)-1i*ueg(:,:,2));
tsubplot(2,2,2); imagesc(gx,gx,log10(uerrg)); showsegment(s);
caxis([-16 0]); colorbar; axis tight; title('log_{10} |error u|');
tsubplot(2,2,3); imagesc(gx,gx,pg); showsegment(s);
caxis([-2 2]); colorbar; axis tight; title('p');
tsubplot(2,2,4); imagesc(gx,gx,log10(abs(pg-peg))); showsegment(s);
caxis([-16 0]); colorbar; axis tight; title('log_{10} error p');
% basic BVP done

% instead use close-evaluation scheme...
ii = s.inside(g.x); g.x = g.x(ii); ug = nan(2*nx^2,1);  % eval only at int pts
tic, ug([ii;ii]) = StoDLP_closeglobal(g,s,mu,tau,'i'); toc  % ii's for 2 cmpts
fprintf('max grid DLP vel cmpt close eval err: %.g\n',max(abs(ug(:)-ueg(:))))
ug = reshape(ug,[nx nx 2]);
figure; set(gcf,'name', 'Sto DLP close u eval');
spg = abs(ug(:,:,1)+1i*ug(:,:,2));   % flow speed on grid
tsubplot(1,2,1); imagesc(gx,gx,spg); showsegment(s);
caxis([0 1.5]); colorbar; axis tight; title('u');
uerrg = abs(ug(:,:,1)+1i*ug(:,:,2)-ueg(:,:,1)-1i*ueg(:,:,2));
tsubplot(1,2,2); imagesc(gx,gx,log10(uerrg)); showsegment(s);
caxis('auto'); colorbar; axis tight; title('log_{10} error u');

% mixed double plus single rep... (not helpful---cond(A) worse---but tests SLP)
A = A + StoSLP(s,s,mu);   % make it the D+S rep (tests S self)
tau = A \ rhs;
fprintf('resid norm %.3g,  density norm %.3g\n',norm(rhs-A*tau),norm(tau))
[upD ppD TpD] = StoDLP(p,s,mu,tau); [upS ppS TpS] = StoSLP(p,s,mu,tau); % eval
up = upD+upS; pp = ppD+ppS; Tp = TpD+TpS;    % vel, pres, trac, at test pts
fprintf('D+S rep: native u error @ test pts:\t\t%.3g\n',max(abs(up-ue(p.x))))
fprintf('\tnative p diff error btw test pts: \t%.3g\n',diff(pp)-diff(pe(p.x)))
poff = mean(pp - pe(p.x));                   % get observed pres offset in rep
Tp = Tp - poff*[-real(p.nx);-imag(p.nx)];      % fix traction's pres offset
fprintf('\tnative traction err @ test pts: \t%.3g\n',max(abs(Tp-Te(p.x,p.nx))))

% again use close-evaluation scheme on same interior pts to test S close...
tic, ug([ii;ii]) = StoDLP_closeglobal(g,s,mu,tau,'i') + StoSLP_closeglobal(g,s,mu,tau,'i'); toc  % ii's for 2 cmpts
fprintf('max grid D+S vel cmpt close eval err: %.g\n',max(abs(ug(:)-ueg(:))))
ug = reshape(ug,[nx nx 2]);
figure; set(gcf,'name', 'Sto D+S close u eval');
spg = abs(ug(:,:,1)+1i*ug(:,:,2));   % flow speed on grid
tsubplot(1,2,1); imagesc(gx,gx,spg); showsegment(s);
caxis([0 1.5]); colorbar; axis tight; title('u');
uerrg = abs(ug(:,:,1)+1i*ug(:,:,2)-ueg(:,:,1)-1i*ueg(:,:,2));
tsubplot(1,2,2); imagesc(gx,gx,log10(uerrg)); showsegment(s);
caxis('auto'); colorbar; axis tight; title('log_{10} error u');
