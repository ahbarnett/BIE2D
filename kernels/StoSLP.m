function [u,p,T] = StoSLP(t,s,mu,dens)
% STOSLP   Evaluate 2D Stokes single-layer velocity, pressure, and traction.
%
% [A,P,T] = StoSLP(t,s,mu) returns dense matrices taking single-layer
%  density values on the nodes of a source curve to velocity, pressure, and
%  traction on the nodes of a target curve. Native quadrature is used, apart
%  from if target=source when A and T are the Nystrom matrices filled using
%  spectrally-accurate quadratures (log-singular for A; smooth diagonal limit
%  for T, ie transpose of DLP).
%
% [u,p,T] = StoSLP(t,s,mu,dens) evaluates the single-layer density dens,
%  returning flow velocity u, pressure p, and target-normal traction T.
%
% The velocity potential is, given fixed viscosity mu>0, and density tau,
% on a source curve gamma, at a target point x,
%
%  u(x) = int_gamma G(x-y) tau(y) ds_y                    x,y in R2
%
%     with 2x2 matrix kernel G(r) = (1/4pi.mu) [ log(1/r).I_2 + r⊗r/|r|^2 ]
%
%  The native quadrature rule of the source curve s is used.
%  For formulae for pressure and traction, see Sec 2.3 of [HW], and [Mar16].
%
% References:
%  [HW]    "Boundary Integral Equations", G. C. Hsiao and W. L. Wendland
%          (Springer, 2008).
%
%  [Mar16] "A fast algorithm for simulating multiphase flows through periodic
%          geometries of arbitrary shape," G. Marple, A. H. Barnett,
%          A. Gillman, and S. Veerapaneni, SIAM J. Sci. Comput. 38(5), B740-772
%          (2016). https://arXiv.org/abs/1510.05616
%
% Inputs: (see setupquad for source & target struct definitions)
%  s = source segment struct with s.x nodes, s.w weights on [0,2pi),
%      s.sp speed function |Z'| at the nodes, and s.tang tangent angles.
%  t = target segment struct with t.x nodes, and t.nx normals if traction needed
%  mu = viscosity
%
% Outputs: (matrix case)
%  A = 2M-by-2N matrix taking density (force vector) to velocity on the
%      target curve. As always for Stokes, ordering is nodes fast,
%      components (1,2) slow, so that A has 4 large blocks A_11, A_12, etc.
%  P = M-by-2N matrix taking density to pressure (scalar) on target nodes.
%  T = 2M-by-2N matrix taking density to normal traction on target nodes.
% Outputs: (density case) all are col vecs (as if N=1).
%
% Notes: 1) Uses whatever self-interaction quadrature Laplace SLP uses
%
% To test call without arguments; also see STOINTDIRBVP
%
% See also: SETUPQUAD, LAPSLP.

% Barnett 6/12/16; T diag limit as in Bowei code SLPmatrixp 6/13/16. 6/27/16.
% vel self-test 4/26/21.
% Todo: * flip args order so mu always first?
if nargin==0, test_StoSLP; return; end  

if numel(mu)~=1, error('mu must be a scalar'); end
if nargout==1
  u = StoSLPmat(t,s,mu);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
  end
elseif nargout==2
  [u p] = StoSLPmat(t,s,mu);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
    p = p * dens;
  end
else
  [u p T] = StoSLPmat(t,s,mu);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
    p = p * dens;
    T = T * dens;
  end
end
%%%%%%

function [A,P,T] = StoSLPmat(t,s,mu)
% Returns native quadrature matrices, or self-evaluation matrices.

self = sameseg(t,s);
N = numel(s.x); M = numel(t.x);
r = bsxfun(@minus, t.x, s.x.');          % C-# displacements mat
irr = 1./(conj(r).*r);                   % 1/r^2, used in all cases below
d1 = real(r); d2 = imag(r);              % worth storing I think
c = 1/(4*pi*mu);              % factor from Hsiao-Wendland book, Ladyzhenskaya
if self
  S = LapSLP(s,s);                    % note includes speed weights
  A = kron((1/2/mu)*eye(2),S);        % prefactor & diagonal log-part blocks
  t1 = real(s.tang); t2 = imag(s.tang); % now add r tensor r part, 4 blocks
  A11 =  d1.^2.*irr; A11(diagind(A11)) = t1.^2;     % diagonal limits
  A12 =  d1.*d2.*irr; A12(diagind(A12)) = t1.*t2;
  A22 =  d2.^2.*irr; A22(diagind(A22)) = t2.^2;
  A = A + bsxfun(@times, [A11 A12; A12 A22], c*[s.w(:)' s.w(:)']); % pref & wei
else                     % distinct src and targ
  logir = -log(abs(r));  % log(1/r) diag block
  A12 = d1.*d2.*irr;     % off diag vel block
  A = c*[logir + d1.^2.*irr, A12;                         % u_x
         A12,                logir + d2.^2.*irr];         % u_y
  A = bsxfun(@times, A, [s.w(:)' s.w(:)']);               % quadr wei
end
if nargout>1           % pressure (no self-eval)
  P = [d1.*irr, d2.*irr];
  P = bsxfun(@times, P, (1/2/pi)*[s.w(:)' s.w(:)']);      % quadr wei
end
if nargout>2           % traction (negative of DLP vel matrix w/ nx,ny swapped)
  rdotn = bsxfun(@times, d1, real(t.nx)) + bsxfun(@times, d2, imag(t.nx));
  rdotnir4 = rdotn.*(irr.*irr); clear rdotn
  A12 = (-1/pi)*d1.*d2.*rdotnir4;
  T = [(-1/pi)*d1.^2.*rdotnir4,   A12;                    % own derivation
     A12,                      (-1/pi)*d2.^2.*rdotnir4];
  if self
    c = -s.cur/2/pi;           % diagonal limit of Laplace DLP
    tx = 1i*s.nx; t1=real(tx); t2=imag(tx);     % tangent vectors on the curve
    T(sub2ind(size(T),1:N,1:N)) = c.*t1.^2;     % overwrite diags of 4 blocks
    T(sub2ind(size(T),1+N:2*N,1:N)) = c.*t1.*t2;
    T(sub2ind(size(T),1:N,1+N:2*N)) = c.*t1.*t2;
    T(sub2ind(size(T),1+N:2*N,1+N:2*N)) = c.*t2.^2;
  end
  T = bsxfun(@times, T, [s.w(:)' s.w(:)']);      % quadr wei
end


%%%%
function test_StoSLP           % only tests vel against reference
b = wobblycurve(1,0.3,5,400);
mu = 0.8;       % viscosity
t.x = 1.5+1i;   % far
densfun = @(t) [1+sin(2+4*t+cos(t));-0.5+cos(3*t)];  % must be 2pi-per, vector
ua = lpevaladapt(t.x,@(x,y,ny) StoSLPvelker(mu,x,y,ny), densfun,b,1e-12);
dens = densfun(b.t);          % eval dens at nodes
u = StoSLP(t,b,mu,dens);
disp(max(abs(u(:)-ua(:))))
