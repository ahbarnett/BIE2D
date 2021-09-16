% test panels (incl close) in 2D for Laplace Green's representation formula.
% Barnett 9/15/21. Helsing-Ojala rule, and potential (no derivs), for now.

clear; a = .3; w = 5;         % smooth wobbly radial shape params
verb = 1;                     % 0 text only; 1 makes a figure
s = wobblycurve(1,a,w,100);   % parametrix descrip of curve (dummy # pts)
Np = 30;                      % # panels (not a convergence study, just a test)
p = 12;                       % panel order
[pa tpan s] = quadr_uniform_panels(s,Np,p);
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
u = @(z) real(fholom(z));                                 % harmonic test func
ux = @(z) real(fpholom(z)); uy = @(z) -imag(fpholom(z));  % partials, sign!
% Since z = x_1 + ix_2 is our packing of coordinates in R2, need to unpack
% to make coord-wise func and grad here...
fprintf('checkgrad: %.3g\n', checkgrad(@(x) u(x(1,:)+1i*x(2,:)), @(x) [ux(x(1,:)+1i*x(2,:)); uy(x(1,:)+1i*x(2,:))]))

% get bdry data u^-, u_n^- for test solution ...
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

% ========= test close interior pt, with i) Helsing-Ojala and ii) sing-swap
s0 = 0.3; dist = 1e-3; t.x = s.Z(s0) - dist*(-1i*s.Zp(s0)/abs(s.Zp(s0)));
if verb, plot(t.x,'r*'); drawnow; end

side = 'i'; closepan = 1.2;    % factor for when to use close

for meth = 'hs'          %  test Helsing-Ojala then singularity-swap
  vt = 0;      % pot at targ
  for i=1:Np   % add in each panel contrib
    a = zpan(i); b = zpan(i+1);     % panel endpts
    if abs(t.x - (a+b)/2) < closepan * abs(b-a)        % near-field of p
    % note: for meth='h', could skip this closepan test. for meth='s' not so.
      D = LapDLP_closepanel(t,pa{i},a,b,side,meth);  % close eval matrix
      S = LapSLP_closepanel(t,pa{i},a,b,side,meth);  % !meth='s' reverts to 'h'
    else
      D = LapDLP(t,pa{i});     % plain rule matrix
      S = LapSLP(t,pa{i});
    end
    jj = (1:p)+p*(i-1);        % indices of pan nodes within global list
    vt = vt + S*unb(jj) - D*ub(jj);      % add contrib of GRF from this pan
  end
  err = vt - u(t.x);      % compare to known u at targ
  fprintf('GRF close int pt val (meth=%s) err = %.3g\n',meth,abs(err))
end
