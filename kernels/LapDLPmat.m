function [A An] = LapDLPmat(t,s)
% [A An] = LapDLPmat(t,s)
% plain double-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg.
% No jump included on self-interaction (ie principal-value integral).
% Self-evaluation for the hypersingular An currently gives inf.

% Barnett 6/12/16 from stuff since 2008.
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
ny = repmat(s.nx.', [M 1]);      % identical rows given by src normals
A = (1/2/pi) * real(ny./d);      % complex form of dipole
if sameseg(t,s)                  % self?
  A(diagind(A)) = -s.cur/4/pi; end           % diagonal term for Laplace
A = A .* repmat(s.w(:)', [M 1]);
if nargout>1     % deriv of double-layer. Not correct for self-interaction.
  csry = conj(ny).*d;              % (cos phi + i sin phi).r
  nx = repmat(t.nx, [1 N]);        % identical cols given by target normals
  csrx = conj(nx).*d;              % (cos th + i sin th).r
  r = abs(d);                      % dist matrix R^{MxN}
  An = -real(csry.*csrx)./(r.^2.^2)/(2*pi);
  An = An .* repmat(s.w(:)', [M 1]);
end

