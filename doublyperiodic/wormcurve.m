function s = wormcurve(a,b,N)
% WORMCURVE   Set up worm-shaped smooth closed curve used in doublyperiodic
%
% s = wormcurve(a,b,N) with a,b shape params, returns segment struct of the
%  form in setupquad. No analytic derivatives used.

% Barnett moved from dpls 6/29/16

s.t = (0:N-1)'/N*2*pi;
s.x = a*cos(s.t)+b*1i*sin(s.t);
s.x = s.x + 0.3i*sin(2*real(s.x));   % spec diff can limit acc in this shape
s = setupquad(s,N);  % uses spectral differentiation, not as accurate as analytic kappa(s)
s.inside = @(z) (real(z)/a).^2+((imag(z)-0.3*sin(2*real(z)))/b).^2 < 1; % yuk
