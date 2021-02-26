function g = LapSLPpotker(x,y,ny)
% LAPSLPPOTKER   kernel func of Laplace single-layer potential, value
%
%   g(x,y) = (1/2pi) log(1/r)   where r:=x-y,  x,y in R2
%
% target x, and source y, are given as complex numbers (C \equiv R2).
%
% Must match LapSLP. Mostly useful as reference, & for adaptive integration.
% ny unused
% Barnett 2/25/21
g = (-1/2/pi) * log(abs(x-y));
