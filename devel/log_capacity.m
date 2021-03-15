% test SLP logarithmic capacity = 1 problems, see if present for Stokes.
% See: Yan-Sloan 1988 JIEA paper, or Hsiao-Wendland book p.11, for Laplace case.
% Barnett 3/15/21

clear
n=50;
mu = 0.7;  % Stokes viscosity param (>0, irrelevant)
fprintf('last three sigma_j(A) for A Laplace or Stokes SLP matrix:\n')

for R = [0.9, 1]    % loop over two circle radii, 2nd being special (log cap=1)
  fprintf('circle R=%g...\n',R)
  s = wobblycurve(R,0,1,n);

  A = LapSLP(s,s);
  ss = svd(A); fprintf('\tLap:'); fprintf('\t%.3g',ss(end-2:end)); fprintf('\n')
  A = A + ones(n,n);    % 1s-matrix fix
  ss = svd(A); fprintf('\tLap+r1:'); fprintf('\t%.3g',ss(end-2:end)); fprintf('\n')

  A = StoSLP(s,s,mu);
  ss = svd(A); fprintf('\tSto:'); fprintf('\t%.3g',ss(end-2:end)); fprintf('\n')
  nv = [real(s.nx);imag(s.nx)];   % col vec of the normals
  %norm(A*nv)                     % check it's R sing vec of A
  A = A + [1;zeros(2*n-1,1)] * nv';    % works
  %A = A + rand(2*n);  % works
  %A = A + [.3*ones(n,1);-0.9*ones(n,1)] * nv';  % fails, orthog to L sing vec
  ss = svd(A); fprintf('\tSto+r1:'); fprintf('\t%.3g',ss(end-2:end)); fprintf('\n')
end
% conclusions:
% * the Laplace log-capacity issue is fixed via rank-1, just like Stokes.
% * there is no log-capacity issue for Stokes SLP, just usual pressure-nullity-1
