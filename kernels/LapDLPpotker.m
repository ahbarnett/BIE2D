function g = LapDLPpotker(x,y,ny)
% LAPDLPPOTKER   kernel func of Laplace double-layer potential, value
%
% g = LapDLPpotker(x,y,ny) evaluates
%
%   g(x;y,ny) = (d/d ny) (1/2pi) log(1/r)   where r:=x-y,  x,y in R2
%
% target x, and source y (normal ny) are given as complex numbers (C \equiv R2).
%
% Mostly useful as reference, & for adaptive integration.

% Barnett 2/25/21
r = abs(x-y);
g = (1/2/pi) * real(conj(ny) .* (x-y)) ./ r.^2;
