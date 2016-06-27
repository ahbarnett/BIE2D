function testStokernels
% TESTSTOKERNELS   plot and test properties of Stokes kernels.
%
% Plot and test 2D Stokes kernels, both velocity and pressure bits,
% check satisfies Stokes PDE, traction is correct, and
% net outflow and force on wall integrated over enclosing circle,

% Barnett 6/27/16 cleaned up from testkernels of 1/29/14-2/11/16.

s.x = 1 + 1.5i;          % src pt, C-#
s.nx = exp(1i*pi/5);     % surface normal (for DLP only), C-#
s.w = 1;                 % src pt quadr wei (dummy for now)
th = pi/6; force = [cos(th);sin(th)];  % src force vector @ angle th
fprintf('src normal (%g,%g), force (%g,%g)\n',...
        real(s.nx),imag(s.nx),force(1),force(2))
fprintf('src normal dot force = %g\n',force(1)*real(s.nx)+force(2)*imag(s.nx))
mu = 0.3;                % viscosity, some random positive value

disp('test (u,p) pairs satisfy Stokes and their traction correct, using finite diff approx at one pt (expect around 9 digits)...')
t.x = 3+2i; t.nx = exp(1i*pi/3);  % target pt and normal for traction eval
fprintf('SLP:'); testPDE(@StoSLPeval,t,s,force,mu)
fprintf('DLP:'); testPDE(@StoDLPeval,t,s,force,mu)

if 1  % vector flow and pressure color plot...
  n = 50; gx = 2*pi*((1:n)-0.5)/n; gy = gx; % target plotting grid
[xx yy] = meshgrid(gx,gy); clear t; t.x = xx(:)+1i*yy(:); clear xx yy
Mt = numel(t.x);

% SLP........ swirling due to pushing a force at src pt, mass is conserved
[u p] = StoSLPeval(t,s,force,mu);
u1=reshape(u(1:n*n),[n n]); u2=reshape(u(n*n+1:end),[n n]); p=reshape(p,[n n]);
figure; imagesc(gx,gy, p);  caxis([-1 1]); hold on; axis xy equal tight;
plot(s.x, 'r*');
sc = 50; quiver(gx,gy, u1,u2, sc./norm([u1(:);u2(:)])); title('Stokes SLP')

% DLP.........  note flow is purely radial, but with angular size variation
figure; %for n = exp(1i*2*pi*(1:100)/100), %v = [real(n);imag(n)]; % anim loop
[u p] = StoDLPeval(t,s,force,mu);
u1=reshape(u(1:n*n),[n n]); u2=reshape(u(n*n+1:end),[n n]); p=reshape(p,[n n]);
imagesc(gx,gy, p); caxis([-1 1]); hold on; axis xy equal tight;
plot(s.x, 'r*');
sc = 200; quiver(gx,gy, u1,u2, sc./norm([u1(:);u2(:)])); title('Stokes DLP')
%hold off; drawnow; end                                           % end anim
end

% Integrate over an enclosing circle...
N=100; R=0.7; t.x = s.x + R*exp(2i*pi*(1:N)'/N); t.nx = exp(2i*pi*(1:N)'/N);
t.w = ones(1,N)*2*pi*R/N;   % targ quadr wei
[u,~,T] = StoSLPeval(t,s,force,mu); Mt = numel(t.x);
u1 = u(1:Mt); u2 = u(Mt+(1:Mt));  % vel
f1 = -T(1:Mt); f2 = -T(Mt+(1:Mt));  % force (-traction dot target nx)
outfl = t.w*(u1.*real(t.nx) + u2.*imag(t.nx));
fprintf('SLP: net outflow %g, \t net force (%g,%g)\n',outfl,t.w*f1,t.w*f2)
[u,~,T] = StoDLPeval(t,s,force,mu); Mt = numel(t.x);
u1 = u(1:Mt); u2 = u(Mt+(1:Mt));  % vel
f1 = -T(1:Mt); f2 = -T(Mt+(1:Mt));  % force (-traction dot target nx)
outfl = t.w*(u1.*real(t.nx) + u2.*imag(t.nx));
fprintf('DLP: net outflow %g, \t net force (%g,%g)\n',outfl,t.w*f1,t.w*f2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function up = velpres(kernel, xtarg, s, dens, mu)
% wrapper to stack vel and pres of a given kernel, at target xtarg. Ie [u1;u2;p]
[u p] = kernel(struct('x',xtarg),s,dens,mu);    % works for Sto{S,D}LPeval
up = [u;p];

function testPDE(kernel,t,s,v,mu)
% uses src s w/ force v and targ t to test SLP or DLP satisfies mu-Stokes PDE,
% and traction correct.
eps = 1e-4;  % finite difference stencil size
rhs = applyStokesPDEs(@(x) velpres(kernel,x,s,v,mu),t.x,mu,eps);
fprintf('\tforce equation err = %g, incompressibility err = %g\n',...
        norm(rhs(1:2)),abs(rhs(3)))
% finite-diff approx to traction @ x,nx...
T = traction(@(x) velpres(kernel,x,s,v,mu),t.x,t.nx,mu,eps);  % send in targ
[~,~,Te] = kernel(t,s,v,mu);   % supposedly exact traction eval
fprintf('\ttraction err = %g\n',norm(T - Te))

function rhs = applyStokesPDEs(f,x,mu,eps)
% check if func f (returning 3 components: u_1, u_2, and p, given C-# loc x)
% obeys mu-Stokes PDE. Outputs the RHS (3 components: f_1, f_2, and rho)
up = nan(3,5);  % three rows are u1,u2,p at each of 5 stencil pts
up(:,1) = f(x); up(:,2) = f(x-eps); up(:,3) = f(x+eps);
up(:,4) = f(x-1i*eps); up(:,5) = f(x+1i*eps);  % do 5-pt stencil evals
gradp = [up(3,3)-up(3,2); up(3,5)-up(3,4)]/(2*eps);
stencil = [-4 1 1 1 1]'/eps^2;
lapu = up(1:2,:)*stencil;
divu = (up(1,3)-up(1,2) + up(2,5)-up(2,4))/(2*eps);
rhs(1:2) = -mu*lapu + gradp;  % Stokes PDE pair
rhs(3)   =  divu;

function T = traction(f,x,nx,mu,eps)
% approx traction vec T of Stokes pair (u_vec, p) given by func f (returning 3
% components: u_x, u_y, and p, given C-# loc x). Barnett 3/2/14
if numel(nx)==1, n=[real(nx); imag(nx)];   % get target normal as 2-by-1
elseif numel(nx)==2, n = nx(:);
  else error('nx must be one C# or 2-component vector'); end
up = nan(3,5);  % three rows are u1,u2,p at each of 5 stencil pts
up(:,1) = f(x); up(:,2) = f(x-eps); up(:,3) = f(x+eps);
up(:,4) = f(x-1i*eps); up(:,5) = f(x+1i*eps); % do evals
p = up(3,1);  % pressure
stress = [up(1:2,3)-up(1:2,2), up(1:2,5)-up(1:2,4)]/(2*eps); % 2-by-2 d_j u_i
stress = stress + stress';         % d_i u_j + d_j u_i
T = -p*n + mu*stress*n;            % std formula
