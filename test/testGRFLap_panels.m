% test panels (incl close) in 2D for Laplace Green's representation formula.
% Barnett 9/14/21. Helsing-Ojala rule, and potential (no derivs), for now.

clear; a = .3; w = 5;         % smooth wobbly radial shape params
verb = 1;
s = wobblycurve(1,a,w,100);   % parametrix descrip of curve (dummy # pts)
Np = 30;             % # panels (not a convergence study, just a test)
p = 12;             % panel order
[stdpan.t, stdpan.w] = gauss(p);
tpan = linspace(0,2*pi,Np+1);   % panel parameter endpoints
s.x = []; s.xp = []; s.w = []; s.nx = []; s.wxp = [];   % override s nodes
for i=1:Np                    % build each panel, oh, also and quadr info in s
  hl = (tpan(i+1)-tpan(i))/2; % param half-size
  t = tpan(i) + hl*(1+stdpan.t);  % global param vals for this pan
  pa{i}.x = s.Z(t); s.x = [s.x; pa{i}.x];
  pa{i}.xp = s.Zp(t); s.xp = [s.xp; pa{i}.xp];
  pa{i}.nx = -1i*pa{i}.xp./abs(pa{i}.xp); s.nx = [s.nx; pa{i}.nx];
  pa{i}.sp = abs(pa{i}.xp);   % speeds (wrt global param)
  pa{i}.w = hl*stdpan.w.'.*pa{i}.sp; s.w = [s.w; pa{i}.w];  % col vec, speed wei
  pa{i}.wxp = hl*stdpan.w.'.*pa{i}.xp; s.wxp = [s.wxp; pa{i}.wxp]; % Helsing wzp
end
zpan = s.Z(tpan);             % panel endpoint locations, in C plane

clear t; t.x = 0.5+0.2i;      % targ: far inside pt (so plain quadr good)
t.nx = exp(1i*pi*1.9);        % arb targ normal

if verb, figure(1); clf; l=0.1; for i=1:Np, plot(pa{i}.x,'k.-'); hold on;
    plot((pa{i}.x + pa{i}.nx*[0 0.1]).','b-');   % show pan nodes, nors, indices
    mx = mean(pa{i}.x); text(real(mx),imag(mx),sprintf('%d',i));
  end, plot(t.x,'g*'); axis equal; drawnow; end

% Let u solve homog PDE (Lap eq) in Omega.
% We test interior GRF formula: u = S u_n^- - D u^-   in Omega, and surface lim

fholom = @(z) z.^4 + exp(z);  fpholom = @(z) 4*z.^3 + exp(z);  % holomorphic
unitDLP = 0;       % 1 overrides with basic tau=-1 const test, no SLP
if unitDLP, fholom = @(z) 1+0*z; fpholom = @(z) 0*z; end
u = @(z) real(fholom(z));                                 % harmonic
ux = @(z) real(fpholom(z)); uy = @(z) -imag(fpholom(z));  % partials, sign!
% Since z = x_1 + ix_2 is our packing of coordinates in R2, need to unpack
% to make coord-wise func and grad here...
fprintf('checkgrad: %.3g\n', checkgrad(@(x) u(x(1,:)+1i*x(2,:)), @(x) [ux(x(1,:)+1i*x(2,:)); uy(x(1,:)+1i*x(2,:))]))

% get interior bdry data u^-, u_n^- ...
ub = u(s.x);
unb = real(s.nx).*ux(s.x) + imag(s.nx).*uy(s.x);    % normal dot grad u

% ========= test far interior point via plain native quadrature eval...
[S DT] = LapSLP(t,s);    % get mat for eval at t, and for directional deriv at t
[D T] = LapDLP(t,s);

% interior GRF formula: u = S u_n^- - D u^-   in Omega
vt = S*unb - D*ub;      % val at targ, by GRF
vnt = DT*unb - T*ub;    % direc-deriv at targ, by GRF
err = vt - u(t.x);      % compare to known u at targ
errn = vnt - (ux(t.x)*real(t.nx) + uy(t.x)*imag(t.nx));   % ..or known deriv
fprintf('GRF far int pt (val,dderiv) err = (%.3g,%.3g)\n',abs(err),abs(errn))

% ========= test close interior pt. with Helsing for now.
t.x = 0.8+0.49i; if verb, plot(t.x,'r*'); drawnow; end   % dist about 5e-3

side = 'i'; closepan = 1.2;    % factor for when to use close
vt = 0;      % pot at targ
for i=1:Np   % add in each panel contrib
  a = zpan(i); b = zpan(i+1);     % panel endpts
  if abs(t.x - (a+b)/2) < closepan * abs(b-a)        % near-field of pan
    D = LapDLP_closepanel(t,pa{i},a,b,side);  % close eval matrix
    S = LapSLP_closepanel(t,pa{i},a,b,side);
  else
    D = LapDLP(t,pa{i});     % plain rule matrix
    S = LapSLP(t,pa{i});
  end
  jj = (1:p)+p*(i-1);        % indices of pan nodes within global list
  vt = vt + S*unb(jj) - D*ub(jj);      % add contrib of GRF from this pan
end
err = vt - u(t.x);      % compare to known u at targ
fprintf('GRF close int pt val err = %.3g\n',abs(err))

