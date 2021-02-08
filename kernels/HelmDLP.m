function [u un] = HelmDLP(k,t,s,dens)
% HELMDLP   Evaluate Helmholtz double-layer potential from curve to targets
%
% This evaluates the 2D Helmholtz double-layer potential for the density tau,
%
%   u(x) = (i/4) int_gamma (d/dn_y) H^{(1)}_0(kr) tau(y) ds_y,  where r:=x-y
%   x,y in R2,
%
%  using the native quadrature rule on the source segment gamma, where point
%  values of tau are given. k is the wavenumber (>0).
%
% [u un] = HelmDLP(k,t,s,dens) evaluates potential and its target-normal
%  derivative.
%
% [A An ...] = HelmDLP(k,t,s) or HelmDLP(k,t,s,[]) returns matrix which maps a
%  density vector to the vector of potentials (A), and possibly target-normal
%  derivatives (An)
%
% If t is the same segment as s, the Kress rule for self-evaluation is used,
%  which assumes s is a global periodic trapezoid rule.
%
% Tested by: TESTGRFLAPHELM

% Crude O(NM) RAM for now. Barnett 2/7/21.
  
if nargout==1
  u = HelmDLPmat(k,t,s);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
  end
elseif nargout==2
  [u un] = HelmDLPmat(k,t,s);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
    un = un * dens;
  end
end
%%%%%%

function [A An] = HelmDLPmat(k,t,s)
% [A An] = HelmDLPmat(k,t,s)
% plain double-layer kernel matrix & targ n-deriv, for Helmholtz in R2.
% k = positive wavenumber, t = target seg (x,nx cols), s = src seg (x,nx cols).
% Kress quadrature used for self-interaction (assumes global quad) in A, An.
% No jump included on self-interaction of value (ie PV integral).
N = numel(s.x);
d = bsxfun(@minus,t.x,s.x.');          % C-# displacements mat
r = abs(d);
cosker = real(bsxfun(@times,conj(s.nx.'),d) ./ r);   % cos ang source nor to d
ss = sameseg(t,s);
if ~ss                                 % plain far-field rule
  A = besselh(1,k*r) .* cosker;
  A = bsxfun(@times, A, (1i*k/4)*s.w(:)');  % prefactor & wei
else
  %***
  A = -log(abs(d)) + circulant(0.5*log(4*sin(pi*(0:N-1)/N).^2));   % peri log
  A(diagind(A)) = -log(s.sp);                       % diagonal limit
  m = 1:N/2-1; Rjn = ifft([0 1./m 2/N 1./m(end:-1:1)])/2;  % Kress Rj(N/2)/4pi
  A = A/N + circulant(Rjn); % includes SLP prefac 1/2pi. Kress peri log matrix L
  A = bsxfun(@times, A, s.sp.');   % do speed factors (2pi/N weights already)
end
if nargout>=2                      % apply T (code from mpspack/@layerpot/T.m)
  if ~ss
    csrx = bsxfun(@times,conj(s.nx.'),d);       % cos targ normals times r
    csry = bsxfun(@times,conj(t.nx),d);       % cos src normals times r
    cc = real(csry).*real(csrx) ./ (r.*r);      % cos phi cos th
    cdor = real(csry.*csrx) ./ (r.*r.*r);   % cos(phi-th) / r
    An = -besselh(1,k*r) .* cdor + k*cc.*besselh(0,k*r);
  else
    % ***  cot self-int?
  end
  An = bsxfun(@times, An, (1i*k/4)*s.w(:)');  % prefactor & wei
end
