% Green's representation formula (GRF) tests for Stokes kernels, vel only.
% Love overdue in this pkg.   Barnett 6/28/21.

clear; a = .3; w = 5;         % smooth wobbly radial shape params
N = 200;                      % convergence param
s = wobblycurve(1,a,w,N);
t.x = 0.4+0.3i;               % targ: far inside pt (so plain quadr good)
t.nx = exp(1i*pi*1.9);        % arb targ normal
mu = 0.7;           % Stokes viscosity of bulk (arb>0)

% Let (u,p) solve homog PDE in Omega; u=vel, p=pressure.
% Test interior GRF formula: u = S T(u,p)^- - D u^-   in Omega, and surface lim
% where T(u,p) is traction.

fprintf('Stokes GRF (mu=%g)...\n',mu)
z0 = 1.5+1.2i;               % source exterior
f0 = [0.6;0.8];              % source strength force vector
u = @(z) StoSLPvelker(mu,z,z0,NaN) * f0;     % Stokeslet at z0, strength f0
T = @(z,nz) StoSLPtracker(mu,z,nz,z0,NaN) * f0;

% get interior bdry data u^-, T^- ...
ub = u(s.x);
Tb = T(s.x,s.nx);

% test interior point via plain native quadrature eval...
[S Spre DT] = StoSLP(t,s,mu);    % get mat for vel, pres, and trac
[D Tpre Tr] = StoDLP(t,s,mu);

% interior GRF formula: u = S T(u,p)^- - D u^-   in Omega
vt = S*Tb - D*ub;      % vel at targ, by GRF
Tt = DT*Tb - Tr*ub;    % trac at targ, by GRF
err = vt - u(t.x);      % compare to known u at targ
errn = Tt - T(t.x,t.nx);      % compare to known u at targ
fprintf('GRF far interior pt (vel,trac) err nrms = (%.3g,%.3g)\n',norm(err),norm(errn))

% Taking interior limit gives    S u_n^- - (D+1/2) u^- = 0  on bdry,
% needs self-eval of the segment...
% GRF resid func, should=0 :
gb = StoSLP(s,s,mu)*Tb - (StoDLP(s,s,mu) + 0.5*eye(2*N))*ub;
fprintf('GRF on-surf vel test rms residual = %.3g\n',norm(gb)/sqrt(N))

% *** could add near-surf target using a close-eval method; todo.
