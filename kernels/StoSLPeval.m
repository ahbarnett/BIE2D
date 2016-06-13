function [u,p,T] = StoSLPeval(t,s,dens,mu)
% STOSLPEVAL  evalutate Stokes single-layer velocity, pressure, and traction.
%
% [u,p,T] = StoSLPeval(t,s,dens,mu) evaluates velocity, pressure, and traction
%  on the nodes of a target curve due to single-layer density values on the
%  nodes of a source curve. Native quadrature is used, unless the target equals
%  the source in which case spectrally-accurate self-evaluation is used.
%  For formulae and normalization see STOSLPMAT.
%
% Inputs: (see setupquad for source & target struct definitions)
%  s = source segment struct with s.x nodes, s.w weights on [0,2pi),
%      s.sp speed function |Z'| at the nodes, and s.tang tangent angles.
%  t = target segment struct with t.x nodes, and t.nx normals if traction needed
%  dens = 2N-by-1 column vector with ordering: nodes fast, components (1,2) slow
%         Can also handle 2N-by-n stack of dens vectors.
%  mu = viscosity
%
% Outputs:
%  u = 2M-by-n velocities on the
%      target curve nodes. As always for Stokes, ordering is nodes fast,
%      components (1,2) slow.
%  p = M-by-n pressure (scalar) on target nodes.
%  T = 2M-by-n normal traction on target nodes.
%
% Notes: 1) Uses whatever self-interaction quadrature Laplace SLP uses
% 2) crude slow for now.
%
% See also: STOSLPMAT.  Tested by: STOINTDIRBVP

% Barnett 6/13/16

if nargout==1
  u = StoSLPmat(t,s,mu) * dens;
elseif nargout==2
  [S P] = StoSLPmat(t,s,mu);
  u = S*dens;
  p = P*dens;
else
  [S P Tr] = StoSLPmat(t,s,mu);
  u = S*dens;
  p = P*dens;
  T = Tr*dens;
end
