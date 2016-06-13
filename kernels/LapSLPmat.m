function [A An] = LapSLPmat(t,s)
% [A An] = LapSLPmat(t,s)
% plain single-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg.
% Kress quadrature used for self-interaction (assumes global quad).
% No jump included on self-interaction of derivative (ie PV integral).

% Barnett 6/12/16 from stuff since 2008.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
if sameseg(t,s)
  A = -log(abs(d)) + circulant(0.5*log(4*sin(pi*(0:N-1)/N).^2));   % peri log
  A(diagind(A)) = -log(s.sp);                       % diagonal limit
  m = 1:N/2-1; Rjn = ifft([0 1./m 2/N 1./m(end:-1:1)])/2;  % Kress Rj(N/2)/4pi
  A = A/N + circulant(Rjn); % includes SLP prefac 1/2pi. Kress peri log matrix L
  A = A .* repmat(s.sp.',[N 1]);  % do speed factors (2pi/N weights already)
else
  A = -(1/2/pi) * log(abs(d)) .* repmat(s.w(:)', [M 1]);
end
if nargout==2                      % apply D^T
  nx = repmat(-t.nx, [1 N]);       % identical cols given by -targ normals
  An = (1/2/pi) * real(nx./d);     % complex form of dipole. Really A1 is An
  if sameseg(t,s)
    An(diagind(An)) = -s.cur/4/pi;  % self? diagonal term for Laplace
  end
  An = An .* repmat(s.w(:)', [M 1]);
end
