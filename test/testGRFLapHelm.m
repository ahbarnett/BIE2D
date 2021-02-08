% Green's representation formula (GRF) tests for Laplace, Helmholtz kernels,
% ie the scalar ones. Love overdue in this pkg.   Barnett 2/7/21.

clear; a = .3; w = 5;         % smooth wobbly radial shape params
N = 200;                      % convergence param
s = wobblycurve(1,a,w,N);
t.x = 0.4+0.3i;               % targ: far inside pt (so plain quadr good)
t.nx = exp(1i*pi*1.9);        % arb targ normal

% Let u solve homog PDE in Omega.
% We test interior GRF formula: u = S u_n^- - D u^-   in Omega, and surface lim

for k=[0 10]                  % ........ wavenumbers to cover Lap, Helm cases

  if k==0
    fprintf('Laplace...\n')
    SLP = @LapSLP;  DLP = @LapDLP;
    fholom = @(z) z.^4 + exp(z);  fpholom = @(z) 4*z.^3 + exp(z);  % holomorphic
    u = @(z) real(fholom(z));                                 % harmonic
    ux = @(z) real(fpholom(z)); uy = @(z) -imag(fpholom(z));  % partials, sign!
  else
    fprintf('Helmholtz (k=%g)...\n',k)
    SLP = @(varargin) HelmSLP(k,varargin{:});       % bake in k
    DLP = @(varargin) HelmDLP(k,varargin{:});
    ang = exp(1i*pi*0.7);                 % arb unit direction (as complex)
    u = @(z) exp(1i*k*real(ang'.*z));     % plane wave
    ux = @(z) 1i*k*real(ang)*u(z); uy = @(z) 1i*k*imag(ang)*u(z);  % partials
  end
  % Since z = x_1 + ix_2 is our packing of coordinates in R2, need to unpack
  % to make coord-wise func and grad here...
  fprintf('checkgrad: %.3g\n', checkgrad(@(x) u(x(1,:)+1i*x(2,:)), @(x) [ux(x(1,:)+1i*x(2,:)); uy(x(1,:)+1i*x(2,:))]))
  
  % get interior bdry data u^-, u_n^- ...
  ub = u(s.x);
  unb = real(s.nx).*ux(s.x) + imag(s.nx).*uy(s.x);    % normal dot grad u
  
  % test interior point via plain native quadrature eval...
  [S DT] = SLP(t,s);    % get mat for eval at t, and for directional deriv at t
  [D T] = DLP(t,s);
  
  % interior GRF formula: u = S u_n^- - D u^-   in Omega
  vt = S*unb - D*ub;      % val at targ, by GRF
  vnt = DT*unb - T*ub;    % direc-deriv at targ, by GRF
  err = vt - u(t.x);      % compare to known u at targ
  errn = vnt - (ux(t.x)*real(t.nx) + uy(t.x)*imag(t.nx));   % ..or known deriv
  fprintf('GRF interior pt (val,dderiv) err = (%.3g,%.3g)\n',abs(err),abs(errn))

  % Taking interior limit gives    S u_n^- - (D+1/2) u^- = 0  on bdry,
  % needs self-eval of the segment...
  gb = SLP(s,s)*unb - (DLP(s,s) + 0.5*eye(N))*ub;   % GRF resid func, should=0
  fprintf('GRF on-surf value test rms residual = %.3g\n',norm(gb)/sqrt(N))

  % If we had self-eval for target-normal deriv of S (ie D^T), or of D
  % (hypersingular T), would also be able to test
  % (D^T-1/2) u_n^- - T u^- = 0  on bdry.
  % But we don't have these yet in kernels
  
  fprintf('\n')
end                           % ............
