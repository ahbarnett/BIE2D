function u = lpevaladapt(x,kerfun,densfun,c,tol)
% LPEVALADAPT  evaluate generic layer potential via adaptive integration
%
% u = lpevaladapt(x,kerfun,densfun,curve,tol)
% Inputs:
%   x - target as complex number
%   kerfun - function handle of form g = kerfun(x,y,ny)
%            for target x, source y,ny.
%   densfun - density function handle, in terms of parameter
%   curve - struct describing a parameterized curve, containing fields:
%           paramdom - (1x2) defines parameter domain of curve
%           Z - (function handle) maps param to complex plane
%           Zp - (function handle) parameter-derivative of Z
%   tol - tolerance (if omitted: 1e-10)
%
% Scalar only for now, and single target
%
% See kernels/LapSLP for a test.

% Barnett 2/25/21. Allow vector densities and potentials 4/26/21.
if nargin<5, tol=1e-10; end
warning('off','MATLAB:integral:MaxIntervalCountReached');

ncomp = numel(densfun(c.paramdom(1)));    % get # vector components

abstol = 1e-12;
% The worst of abstol and |Q|*reltol is taken; according to "help integral".
% If push to 1e-13 or below, it requests >1e5 nodes, v slow :(

% integrand = kernel(x,y,ny) * speed (Jacobean from param to arclength) * dens:
if ncomp==1
  % scalar case (t can be array of ordinates)
  u = integral(@(t) kerfun(x,c.Z(t),-1i*c.Zp(t)./abs(c.Zp(t))).*abs(c.Zp(t)).*densfun(t),c.paramdom(1),c.paramdom(2),'abstol',abstol,'reltol',tol);
else
  % vector case (t now is only ever one number)
  u = integral(@(t) kerfun(x,c.Z(t),-1i*c.Zp(t)/abs(c.Zp(t)))*abs(c.Zp(t))*densfun(t),c.paramdom(1),c.paramdom(2),'abstol',abstol,'reltol',tol,'ArrayValued',1);
end
