% Basic solver for Stokes exterior velocity (Dirichlet; no-slip) BVP on smooth
% curve, via completed D+S formulation, and using iterative solution.
% No known solution; convergence study.
% Note: by "Stokes paradox" the velocity induced by the obstacle grows as
% log(r) at infinity; so this is non-physical.
%
% Barnett 7/7/16, modifed from StointDirBVP.m
%
% References: see Sec 2.3 of
%  [HW]    "Boundary Integral Equations", G. C. Hsiao and W. L. Wendland
%          (Springer, 2008).

clear;
mu = 0.7;              % viscosity
a = .3; w = 5;         % smooth wobbly radial shape params

p = []; p.x = [2+1.5i];                     % 1 test pt

uinc = @(x) [1+0*x;0*x]; pinc = @(x) 0*x;   % incident uniform flow

Ns = 40:20:200; ut = nan(2,numel(Ns));
for i=1:numel(Ns), N = Ns(i); % ----------- N-convergence
  s = wobblycurve(a,w,N); s.a = mean(s.x); %figure; showsegment({s p})
  A = eye(2*N)/2 + StoDLP(s,s,mu) + StoSLP(s,s,mu);   % Nystrom, ext JR
  rhs = -uinc(s.x);                          % no-slip velocity data
  %tau = A \ rhs;                          % direct, or...
  tau = gmres(A,rhs,[],1e-14,size(A,1));   % ...iterative
  ut(:,i) = StoDLP(p,s,mu,tau) + StoSLP(p,s,mu,tau);
  fprintf('\tN=%d\tu @ test pt =(%.16g,%.16g)\n',N,ut(1,i),ut(2,i))
end
errs = sqrt(sum((ut-repmat(ut(:,end),[1 numel(Ns)])).^2,1)); % 2-norm u errs
figure; semilogy(Ns,errs,'+-'); xlabel('N'); ylabel('2-norm of vel err');
title('ext no-slip BVP: pointwise far-field N-convergence'); drawnow
%print -dpng stoextnoslip_wobbly_w8_a.5_Nconv.png

% plot solution on grid & soln errors...
nx = 50; gx = 3*linspace(-1,1,nx); gy=gx; ng = nx^2;
[xx yy] = meshgrid(gx); zz = xx(:)+1i*yy(:); clear xx yy
ii = ~s.inside(zz); g.x = zz(ii);
ug = nan(2*ng,1); pg = nan(ng,1);
tic
[uD pD] =  StoDLP_closeglobal(g,s,mu,tau,'e');
[uS pS] =  StoSLP_closeglobal(g,s,mu,tau,'e');
toc
ug([ii;ii]) = uinc(g.x) + uD+uS; pg(ii) = pinc(g.x) + pD+pS;  % total field
u1 = reshape(ug(1:ng),[nx nx]); u2 = reshape(ug(ng+1:end),[nx nx]);
pg = reshape(pg,[nx nx]);
figure; set(gcf,'name', 'ext no-slip BVP total u,p (closeglobal eval)');
contourf(gx,gy,pg,[min(pg(:)):.5:max(pg(:))]);
colormap(jet(256)); axis equal tight; hold on; plot([s.x;s.x(1)],'k-');
quiver(real(zz(ii)),imag(zz(ii)),u1(ii),u2(ii),2.0,'k-');  % show vec field

%print -dpng stoextnoslip_wobbly_w8_a.5_soln.png