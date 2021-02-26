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

% Barnett 2/25/21
if nargin<5, tol=1e-10; end

% integrand = kernel(x,y,ny) * speed (Jacobean from param to arclength) * dens
u = integral(@(t) kerfun(x,c.Z(t),-1i*c.Zp(t)./abs(c.Zp(t))).*abs(c.Zp(t)).*densfun(t),c.paramdom(1),c.paramdom(2),'abstol',tol,'reltol',tol);
