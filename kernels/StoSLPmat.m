function [A,T] = StoSLPmat(t,s,mu)
% STOSLPMAT   Matrix for Stokes single-layer vel potential and its traction.
%
% Returns 2M-by-2N matrix A from source force vector density to flow vector
%  components at the target, and optionally the same-sized matrix T evaluating
%  traction at the targets. Native quadrature is used.
%
% The normalization is as in Ladyzhenskaya or Hsiao-Wendland book, i.e. a
%  prefactor of 1/(4.pi.mu).
%
% Inputs: (see setupquad for source & target struct definitions)
%  s = source segment struct with s.x nodes, s.w weights on [0,2pi),
%      s.sp speed function |Z'| at the nodes, and s.tang tangent angles.
%  t = target segment struct with t.x nodes, and t.nx normals if traction needed
%  mu = viscosity
%
% Uses whatever self-interaction quadrature Laplace SLP uses
%
% See also: LAPSLPMAT

% Barnett 6/12/16
N = numel(s.x); M = numel(t.x);
r = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
irr = 1./(conj(r).*r);         % 1/r^2, used in all cases below
d1 = real(r); d2 = imag(r);
c = 1/(4*pi*mu);              % factor from Hsiao-Wendland book, Ladyzhenskaya
if sameseg(t,s)
  S = LapSLPmat(s,s);                   % note includes speed weights.
  A = (1/2/mu) * kron(eye(2),S);        % prefactor & diagonal log-part blocks
  t1 = real(s.tang); t2 = imag(s.tang);  % now do r tensor r part...
  A11 =  d1.^2.*irr; A11(diagind(A11)) = t1.^2;     % diagonal limits
  A12 =  d1.*d2.*irr; A12(diagind(A12)) = t1.*t2;
  A22 =  d2.^2.*irr; A22(diagind(A22)) = t2.^2;
  A = A + c*[A11 A12; A12 A22].*repmat([s.w(:)' s.w(:)'], [2*M 1]); % pref & wei
else                     % distinct src and tar
  logir = -log(abs(r));  % log(1/r) diag block
  A12 = d1.*d2.*irr;     % off diag vel block
  A = c*[logir + d1.^2.*irr, A12;                         % u_x
         A12,                logir + d2.^2.*irr];         % u_y
  A = A .* repmat([s.w(:)' s.w(:)'], [2*M 1]);            % quadr wei
end
if nargout>1           % traction (negative of DLP vel matrix w/ nx,ny swapped)
  rdotn = d1.*repmat(real(t.nx), [1 N]) + d2.*repmat(imag(t.nx), [1 N]);
  rdotnir4 = rdotn.*(irr.*irr); clear rdotn
  A12 = -(1/pi)*d1.*d2.*rdotnir4;
  T = [-(1/pi)*d1.^2.*rdotnir4,   A12;                    % own derivation
     A12,                      -(1/pi)*d2.^2.*rdotnir4];
  T = T .* repmat([s.w(:)' s.w(:)'], [2*M 1]);            % quadr wei
end
