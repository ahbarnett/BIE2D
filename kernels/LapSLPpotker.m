function g = LapSLPpotker(x,y,ny)
% LAPSLPPOTKER   kernel func of Laplace single-layer potential, value
%
% g = LapSLPpotker(x,y,ny) evaluates
%
%   g(x,y) = (1/2pi) log(1/r)   where r:=x-y,  x,y in R2
%
% target x, and source y, are given as complex numbers (C \equiv R2).
% ny is unused.
%
% Mostly useful as reference, & for adaptive integration

% Barnett 2/25/21
g = (-1/2/pi) * log(abs(x-y));
