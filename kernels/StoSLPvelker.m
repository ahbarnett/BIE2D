function G = StoSLPvelker(mu,x,y,ny)
% STOSLPVELKER   kernel func of Stokes single-layer potential, velocity
%
% G = StoSLPvelker(mu,x,y,ny) evaluates
%
%   G(x,y) = (1/4pi.mu) [ log(1/r).I_2 + r⊗r/|r|^2 ],  where r:=x-y,  x,y in R2
%
% target x, and source y, are given as complex numbers (C \equiv R2).
% ny is unused.
%
% Mostly useful as reference, & for adaptive integration

% Barnett 4/26/21
r = x-y;
irr = 1./(conj(r).*r);
l = 0.5*log(irr);
G12 = real(r).*imag(r).*irr;
G = (1/4/pi/mu) * [l + real(r).*real(r).*irr,   G12;
                   G12,                         l + imag(r).*imag(r).*irr];
