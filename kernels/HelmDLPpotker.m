function g = HelmDLPpotker(k,x,y,ny)
% HELMDLPPOTKER   kernel func of Helmholtz double-layer potential, value
%
% g = HelmDLPpotker(k,x,y,ny) evaluates DLP kernel for wavenumber-k Helmholtz,
%
%   g_k(x;y,ny) = (d/d ny) (i/4) H^(1)_0(kr)   where r:=x-y,  x,y in R2
%
% target x, and source y, and unit normal are given as complex numbers (C \equiv R2).
%
% Mostly useful as reference, & for adaptive integration.

% Barnett 2/25/21

r = abs(x-y);
costh = real(conj(ny).*(x-y)) ./ r;
g = (1i*k/4) * costh .* besselh(1,k*r);
