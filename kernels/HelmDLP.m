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
  A = besselh(1,k*r) .* cosker;
  A = bsxfun(@times, A, (1i*k/4)*s.w(:)');  % prefactor & wei
  A(diagind(A)) = -(s.cur.*s.w)/(4*pi);    % should be O(1/N^3)
  %***   wrong for now
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
