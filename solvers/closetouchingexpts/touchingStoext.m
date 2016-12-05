% Exploring Stokes close-to-touching curve GMRES convergence, two noslip discs.
% A = dense matrix from completed D+S formulation.
% C = blocks of A involving "close" region only, unit diag elsewhere, A = B+C.
% Limited by use of global quadr for now, which is slow for 
%
% Note: by "Stokes paradox" the velocity induced by the obstacle grows as
% log(r) at infinity; so this is non-physical, but not the issue.
%
% Barnett 12/3/16, modifed from StoextDirBVPconv.m

clear;
mu = 0.7;              % viscosity
ep = 1e-6;             % separation of discs
closequad = 1;         % 0 - plain Nystrom, 1 - spectral close quadr for A12 etc
          % (note 1 is slow eg 10s @ N=300, stupid O(N^3) close matrix fill)
precond = 1;           % 0 - GMRES for whole sys, 1 - direct precond of close
tol = 1e-12;           % GMRES rel resid tol

p = []; p.x = [2+1.5i];                     % a distant test pt for soln conv

uinc = @(x) [0.5+0*x;1+0*x]; pinc = @(x) 0*x;   % incident uniform flow

Ns = 50:50:300;        % list of nodes per disc in convergence study
ut = nan(2,numel(Ns));
for i=1:numel(Ns), N = Ns(i); % ----------- N-convergence
  s1 = wobblycurve(1-ep/2,0,1,N); s1.x=s1.x+1; s1.inside = @(z) s1.inside(z-1); s1.a = mean(s1.x);
  s2 = wobblycurve(1-ep/2,0,1,N); s2.x=s2.x-1; s2.inside = @(z) s2.inside(z+1); s2.a = mean(s2.x);
  %figure; showsegment({s1 s2 p}), stop
  % fill four blocks btw 2 obstacles...
  A11 = eye(2*N)/2 + StoDLP(s1,s1,mu) + StoSLP(s1,s1,mu);
  A22 = eye(2*N)/2 + StoDLP(s2,s2,mu) + StoSLP(s2,s2,mu);
  if closequad
    A21 = StoDLP_closeglobal(s2,s1,mu,[],'e') + StoSLP_closeglobal(s2,s1,mu,[],'e');
    A12 = StoDLP_closeglobal(s1,s2,mu,[],'e') + StoSLP_closeglobal(s1,s2,mu,[],'e');
  else
    A21 = StoDLP(s2,s1,mu) + StoSLP(s2,s1,mu);
    A12 = StoDLP(s1,s2,mu) + StoSLP(s1,s2,mu);
  end
  A = [A11 A12; A21 A22];
  rhs = -[uinc(s1.x); uinc(s2.x)];          % no-slip velocity data
  if precond
    d = 0.05;                % x coord cutoff for 'close' region
    J1 = find(real(s1.x)<d); J1 = [J1;J1+N];   % 'close' indices on s1
    J2 = find(real(s2.x)>-d); J2 = [J2;J2+N];
    J = [J1; J2+2*N];  % index set within full lin sys
    C = eye(4*N); iC=C; C(J,J) = A(J,J); B = A-C; iC(J,J) = inv(C(J,J));
    %norm(C*iC-eye(4*N)),     norm(iC*C-eye(4*N))        % check inverse
    y = gmres(@(x) x + B*(iC*x),rhs,[],tol,4*N);
    %y = (eye(4*N)+B*iC) \ rhs;                          % direct, or...  
    tau = iC*y;
  else
    %tau = A \ rhs;                          % direct, or...  
    tau = gmres(A,rhs,[],tol,size(A,1));   % ...iterative
  end
  tau1 = tau(1:2*N); tau2 = tau(2*N+1:end);  % pull out each curve tau
  ut(:,i) = StoDLP(p,s1,mu,tau1) + StoSLP(p,s1,mu,tau1) + StoDLP(p,s2,mu,tau2) + StoSLP(p,s2,mu,tau2);
  fprintf('\tN=%d\tu @ test pt =(%.16g,%.16g)\n',N,ut(1,i),ut(2,i))
end                    % ---------------
errs = sqrt(sum((ut-repmat(ut(:,end),[1 numel(Ns)])).^2,1)) % 2-norm u errs
%figure; semilogy(Ns,errs,'+-'); xlabel('N'); ylabel('2-norm of vel err');
%title('ext no-slip BVP: pointwise far-field N-convergence'); drawnow
%print -dpng stoextnoslip_wobbly_w8_a.5_Nconv.png

figure; plot(tau,'.-'); title('density');

figure; l=eig(A); plot(real(l),imag(l),'+'); axis equal; title('spec(A)');
if precond, hold on; l=eig(eye(4*N)+B*iC); plot(real(l),imag(l),'r+'); title('spec(A) in blue, spec(I-BC^{-1}) in red'); end

if 1  % plot solution on grid & soln errors...
nx = 80; gx = 3*linspace(-1,1,nx); gy=gx; ng = nx^2;
[xx yy] = meshgrid(gx); zz = xx(:)+1i*yy(:); clear xx yy
ii = ~s1.inside(zz) & ~s2.inside(zz); % outside both
g.x = zz(ii);
ug = nan(2*ng,1); pg = nan(ng,1);
tic
[uD1 pD1] =  StoDLP_closeglobal(g,s1,mu,tau1,'e');
[uD2 pD2] =  StoDLP_closeglobal(g,s2,mu,tau2,'e');
[uS1 pS1] =  StoSLP_closeglobal(g,s1,mu,tau1,'e');
[uS2 pS2] =  StoSLP_closeglobal(g,s2,mu,tau2,'e');
toc
ug([ii;ii]) = uinc(g.x) + uD1+uS1 + uD2+uS2; pg(ii) = pinc(g.x) + pD1+pS1+pD2+pS2;  % total field
u1 = reshape(ug(1:ng),[nx nx]); u2 = reshape(ug(ng+1:end),[nx nx]); % cpts
pg = reshape(pg,[nx nx]);
figure; set(gcf,'name', 'ext no-slip BVP total u,p (closeglobal eval)');
contourf(gx,gy,pg,[min(pg(:)):0.5:max(pg(:))]);
colormap(jet(256)); axis equal tight; hold on; plot([s1.x;s1.x(1)],'k-');plot([s2.x;s2.x(1)],'k-');
quiver(real(zz(ii)),imag(zz(ii)),u1(ii),u2(ii),2.0,'k-');  % show vec field
title('arrows = u vel, contours = p pres')
end
