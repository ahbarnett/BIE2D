% Basic solver for Laplace interior Dirichlet BVP on smooth curve.
% Also serves to test the native SLP and DLP matrix and eval routines.
% Barnett 6/12/16
%
% References: see Ch. 6 of this:
%
% [LIE]  R. Kress, "Linear Integral Equations", vol. 82 of Appl. Math. Sci.,
%        Springer, second ed., 1999

a = .3; w = 5;         % smooth wobbly radial shape params
reps = 'dm';   % what expts to do
N = 200;
s = wobblycurve(a,w,N);

p = []; p.x = .2+.3i; p.nx = exp(1.9i);   % test pt including a direction
fholom = @(z) exp(1i*(z+1));        % holomorphic exact interior soln
fpholom = @(z) 1i*exp(1i*(z+1));    % its complex derivative

f = @(z) real(fholom(z)); imf = @(z) imag(fholom(z)); % real part (true f), etc
fx = @(z) real(fpholom(z)); fy = @(z) -imag(fpholom(z)); % note sign!

A = -eye(N)/2 + LapDLPmat(s,s);       % discretized integral operator

*** to finish - just do far field test pt,
and do mixed D+S rep to test S.
