% Basic solver for Laplace interior Dirichlet BVP on smooth curve.
% Also serves to test the native SLP and DLP matrix and (close) eval routines.
% Barnett 6/12/16-6/27/16
%
% References: see Ch. 6 of this:
%
% [LIE]  R. Kress, "Linear Integral Equations", vol. 82 of Appl. Math. Sci.,
%        Springer, second ed., 1999

clear; a = .3; w = 5;         % smooth wobbly radial shape params
N = 200; s = wobblycurve(a,w,N);

p = []; p.x = .2+.3i; p.nx = exp(1.9i);   % test pt including a direction
figure; showsegment({s p})
fholom = @(z) exp(1i*(z+1));        % holomorphic exact interior soln
fpholom = @(z) 1i*exp(1i*(z+1));    % its complex derivative

f = @(z) real(fholom(z));           % use real part as known Laplace soln
fx = @(z) real(fpholom(z)); fy = @(z) -imag(fpholom(z)); % partials, note sign!


% ------- double-layer representation...
A = -eye(N)/2 + LapDLP(s,s);       % Nystrom discretized integral operator
rhs = f(s.x);
tau = A \ rhs;
fprintf('Lap int Dir BVP resid norm %.3g,  density norm %.3g\n',norm(rhs-A*tau),norm(tau))
[up unp] = LapDLP(p,s,tau);        % potential and targ n-deriv at test pt
fnp = fx(p.x)*real(p.nx) + fy(p.x)*imag(p.nx);   % exact targ n-deriv
fprintf('D rep: native u and un errors @ test pt: \t%.3g\t%.3g \n',up-f(p.x),unp-fnp)

% plot solution on grid & soln errors...
nx = 200; gx = max(abs(s.x))*linspace(-1,1,nx);
[xx yy] = meshgrid(gx); g.x = xx(:)+1i*yy(:);
tic; ug = LapDLP(g,s,tau); toc          % u on grid, expensive bit for now
fg = f(g.x);                    % known soln on grid
ug = reshape(ug,[nx nx]); fg = reshape(fg,[nx nx]);  % shape arrays for plot
figure; set(gcf,'name', 'Lap DLP native eval');
tsubplot(1,2,1); imagesc(gx,gx,ug); showsegment(s);
caxis([-1 2]); colorbar; axis tight; title('u');
tsubplot(1,2,2); imagesc(gx,gx,log10(abs(ug-fg))); showsegment(s);
caxis([-16 0]); colorbar; axis tight; title('log_{10} error u');
% basic BVP done ------

% instead use close-evaluation scheme...
ii = s.inside(g.x); g.x = g.x(ii); ug = nan*ug;  % eval only at interior pts
tic; ug(ii) = LapDLP_closeglobal(g,s,tau,'i'); toc
fprintf('max grid DLP pot close eval err:\t%.3g\n',max(abs(ug(:)-fg(:))))
figure; set(gcf,'name', 'Lap DLP close eval');
tsubplot(1,2,1); imagesc(gx,gx,ug); showsegment(s);
caxis([-1 2]); colorbar; axis tight; title('u');
tsubplot(1,2,2); imagesc(gx,gx,log10(abs(ug-fg))); showsegment(s);
caxis('auto'); colorbar; axis tight; title('log_{10} error u');


% --------- Mixed double plus single rep...
% (not helpful---cond(A) worse---but it tests SLP)
A = A + LapSLP(s,s);   % make it the D+S rep
tau = A \ rhs;
fprintf('resid norm %.3g,  density norm %.3g\n',norm(rhs-A*tau),norm(tau))
[up unp] = LapDLP(p,s,tau); [vp vnp] = LapSLP(p,s,tau);
up = up+vp; unp = unp+vnp;
fprintf('D+S rep: native u and un errors @ test pt: \t%.3g\t%.3g \n',up-f(p.x),unp-fnp)

% again use close-evaluation scheme on same interior pts only, now D+S...
tic; ug(ii) = LapDLP_closeglobal(g,s,tau,'i') + LapSLP_closeglobal(g,s,tau,'i'); toc
fprintf('max grid D+S pot close eval err:\t%.3g\n',max(abs(ug(:)-fg(:))))
figure; set(gcf,'name', 'Lap D+S close eval');
tsubplot(1,2,1); imagesc(gx,gx,ug); showsegment(s);
caxis([-1 2]); colorbar; axis tight; title('u');
tsubplot(1,2,2); imagesc(gx,gx,log10(abs(ug-fg))); showsegment(s);
caxis('auto'); colorbar; axis tight; title('log_{10} error u');
% done

