function [t0 c cp] = panel_preimage(z,t,z0)
% PANEL_PREIMAGE   solve for preimage of complex point under panel map.
%
% t0 = panel_preimage(z,t,z0) returns t0 preimage such that z0 = Z(t0) where
%  Z : C -> C is the analytic map that took the set of standard nodes t to
%  the set of nodes z.
%
% [t0 c cp] = panel_preimage(z,t,z0) also returns the coeffs (in polyval order)
%  of the approximation to the map Z and its complex derivative respectively.
%  c and cp do not depend on z0.

% to do: vectorize over z0, t0

if nargin==0, test_panel_preimage; return; end

% solve monomial rep of map via Vandermonde (indep of z0)
t = t(:); z = z(:); p = numel(t);
V = ones(p,p);
for k=2:p
  V(:,k) = t.*V(:,k-1);
end
c = V\z;                   % monomial coeffs
cp = c(2:end).*(1:p-1)';   % monomial coeffs of deriv of map
c=flipud(c); cp=flipud(cp);   % order max power to const term, for polyval

% Newton to solve for t0 in Z(t0)-z0 = 0.  todo: check vectorizes over z0, t0
maxit = 20;
zcen = (z(p)+z(1))/2; zsc = (z(p)-z(1))/(t(p)-t(1));
t0 = (z0-zcen)/zsc;    % initial guess
for i=1:maxit
  t0old = t0;
  t0 = t0 - (polyval(c,t0) - z0) ./ polyval(cp,t0);
  if max(abs(t0-t0old)) < 1e-15             % not rel since on [-1,1] scale
    break;
  end
end
%%%%%%%%

function test_panel_preimage
Z = @(t) exp((0.2+0.9i)*t);   % analytic but not an exact poly
%Z = @(t) 2 + 1i*t + 0.1*t.^2;  % exact poly
p=16;
t = gauss(p);
z = Z(t);
t0 = 0.7+0.6i;
z0 = Z(t0);
abs(t0 - panel_preimage(z,t,z0))
%figure; plot(z,'k.'); hold on; plot(z0,'r*'); axis equal;
