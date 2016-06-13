function s = wormcurve(a,b,N)
% WORMCURVE   Return smooth curve quadrature for testing periodic solvers
%
% s = wormcurve(a,b,N) where a,b are shape params, and N is the number of nodes.
%  The struct s also includes s.inside which is a function handle returning
%  true when its argument is inside the curve.

s.t = (1:N)'/N*2*pi; s.x = a*cos(s.t)+b*1i*sin(s.t);
s.x = s.x + 0.3i*sin(2*real(s.x));   % spec diff can limit acc in this shape
s.inside = @(z) (real(z)/a).^2+((imag(z)-0.3*sin(2*real(z)))/b).^2 < 1; % yuk
% uses spectral differentiation, not as acc as analytic kappa(s)...
s = setupquad(s);
