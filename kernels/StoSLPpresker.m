function Gp = StoSLPpresker(mu,x,y,ny)
% STOSLPPRESKER   kernel func of Stokes single-layer potential, pressure
%
% Gp = StoSLPpresker(mu,x,y,ny) evaluates
%
%   G_p(x,y) = (1/2pi) r/|r|^2,     where  r:=x-y,  x,y in R2
%
% target x, and source y, are given as complex numbers (C \equiv R2).
% ny and mu are unused.
%
% Mostly useful as reference, & for adaptive integration

% Barnett 6/14/21
r = x-y;
irr = 1./(conj(r).*r);
Gp = (1/2/pi) * [real(r).*irr, imag(r).*irr];
