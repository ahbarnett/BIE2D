function [u,p,T] = StoDLP(t,s,mu,dens)
% STODLP   Evaluate 2D Stokes single-layer velocity, pressure, and traction.
%
% [A,P,T] = StoDLP(t,s,mu) returns dense matrices taking double-layer
%  density values on the nodes of a source curve to velocity, pressure, and
%  traction on the nodes of a target curve. Native quadrature is used, apart
%  from when target=source when A is the self-interaction Nystrom matrix using
%  the smooth diagonal-limit formula for DLP.
%
% [u,p,T] = StoDLP(t,s,mu,dens) evaluates the double-layer density dens,
%  returning flow velocity u, pressure p, and target-normal traction T.
%
% The velocity potential is, given fixed viscosity mu>0, and density tau,
% on a source curve gamma, at a target point x,
%
%  u(x) = int_gamma G(x-y) tau(y) ds_y                    x,y in R2
%
%     with 2x2 matrix kernel G(r) = (1/pi) (r⊗r)(r.n_y)/|r|^4
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
%
% Notes: 1) no self-evaluation for P or T.
%
% To test call without arguments; also see STOINTDIRBVP
%
% See also: SETUPQUAD

% Barnett 6/13/16, 6/27/16
% vel self-test 4/26/21.
% Todo: * flip args order so mu always first?
if nargin==0, test_StoDLP; return; end  

if numel(mu)~=1, error('mu must be a scalar'); end
if nargout==1
  u = StoDLPmat(t,s,mu);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
  end
elseif nargout==2
  [u p] = StoDLPmat(t,s,mu);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
    p = p * dens;
  end
else
  [u p T] = StoDLPmat(t,s,mu);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
    p = p * dens;
    T = T * dens;
  end
end
%%%%%%

function [A,P,T] = StoDLPmat(t,s,mu)
% Returns native quadrature matrices, or self-evaluation matrices.
% some repmat->bsxfun 6/28/16

N = numel(s.x); M = numel(t.x);
r = bsxfun(@minus,t.x,s.x.');                 % C-# displacements mat
irr = 1./(conj(r).*r);    % 1/r^2, used in all cases below
d1 = real(r); d2 = imag(r);
%rdotny = d1.*repmat(real(s.nx)', [M 1]) + d2.*repmat(imag(s.nx)', [M 1]);
rdotny = bsxfun(@times, d1, real(s.nx)') + bsxfun(@times, d2, imag(s.nx)');
rdotnir4 = rdotny.*(irr.*irr); if nargout==1, clear rdotny; end
A12 = (1/pi)*d1.*d2.*rdotnir4;  % off diag vel block
A = [(1/pi)*d1.^2.*rdotnir4,   A12;                     % Ladyzhenskaya
     A12,                      (1/pi)*d2.^2.*rdotnir4];
if sameseg(t,s)         % self-int DLP velocity diagonal correction
  c = -s.cur/2/pi;           % diagonal limit of Laplace DLP
  tx = 1i*s.nx; t1=real(tx); t2=imag(tx);     % tangent vectors on the curve
  A(sub2ind(size(A),1:N,1:N)) = c.*t1.^2;     % overwrite diags of 4 blocks
  A(sub2ind(size(A),1+N:2*N,1:N)) = c.*t1.*t2;
  A(sub2ind(size(A),1:N,1+N:2*N)) = c.*t1.*t2;
  A(sub2ind(size(A),1+N:2*N,1+N:2*N)) = c.*t2.^2;
end
A = bsxfun(@times, A, [s.w(:)' s.w(:)']);            % quadr wei
if nargout>1           % pressure of DLP (no self-int)
%  P = [(mu/pi)*(-repmat(real(s.nx)', [M 1]).*irr + 2*rdotnir4.*d1), ...
%       (mu/pi)*(-repmat(imag(s.nx)', [M 1]).*irr + 2*rdotnir4.*d2)  ];
  P = [bsxfun(@times, real(-s.nx)',irr) + 2*rdotnir4.*d1, ...
       bsxfun(@times, imag(-s.nx)',irr) + 2*rdotnir4.*d2  ];
  P = bsxfun(@times, P, (mu/pi)*[s.w(:)' s.w(:)']);            % quadr wei
end
if nargout>2           % traction, my derivation of formula (no self-int)
  nx1 = repmat(real(t.nx), [1 N]); nx2 = repmat(imag(t.nx), [1 N]);
  rdotnx = d1.*nx1 + d2.*nx2;
  ny1 = repmat(real(s.nx)', [M 1]); ny2 = repmat(imag(s.nx)', [M 1]);
  dx = rdotnx.*irr; dy = rdotny.*irr; dxdy = dx.*dy;
  R12 = d1.*d2.*irr; R = [d1.^2.*irr, R12; R12, d2.^2.*irr];
  nydotnx = nx1.*ny1 + nx2.*ny2;
  T = R.*kron(ones(2), nydotnx.*irr - 8*dxdy) + kron(eye(2), dxdy);
  T = T + [nx1.*ny1.*irr, nx1.*ny2.*irr; nx2.*ny1.*irr, nx2.*ny2.*irr] + ...
      kron(ones(2),dx.*irr) .* [ny1.*d1, ny1.*d2; ny2.*d1, ny2.*d2] + ...
      kron(ones(2),dy.*irr) .* [d1.*nx1, d1.*nx2; d2.*nx1, d2.*nx2];
  T = bsxfun(@times, T, (mu/pi)*[s.w(:)' s.w(:)']);    % prefac & quadr wei
end


%%%%
function test_StoDLP           % only tests vel against reference
b = wobblycurve(1,0.3,5,400);
mu = 0.8;       % viscosity
t.x = 1.5+1i;   % far
densfun = @(t) [1+sin(2+4*t+cos(t));-0.5+cos(3*t)];  % must be 2pi-per, vector
ua = lpevaladapt(t.x,@(x,y,ny) StoDLPvelker(mu,x,y,ny), densfun,b,1e-12);
dens = densfun(b.t);          % eval dens at nodes
u = StoDLP(t,b,mu,dens);
disp(max(abs(u(:)-ua(:))))
