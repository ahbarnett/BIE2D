function g = HelmSLPpotker(k,x,y,ny)
% HELMSLPPOTKER   kernel func of Helmholtz single-layer potential, value
%
% g = HelmSLPpotker(k,x,y,ny) evaluates wavenumber-k Helmholtz fundamental soln
%
%   g_k(x,y) = (i/4) H^(1)_0(kr)   where r:=x-y,  x,y in R2
%
% target x, and source y, are given as complex numbers (C \equiv R2).
% ny unused.
%
% Mostly useful as reference, & for adaptive integration.

% Barnett 2/25/21
g = (1i/4) * besselh(0,k*abs(x-y));
