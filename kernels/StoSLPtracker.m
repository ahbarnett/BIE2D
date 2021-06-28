function G = StoSLPtracker(mu,x,nx,y,ny)
% STOSLPTRACKER   kernel func of Stokes single-layer potential, traction
%
% G = StoSLPtracker(mu,x,nx,y,ny) evaluates
%
%   G(x,y) = (-1/pi) (r⊗r) (r.n_y)/ |r|^4 ,  where r:=x-y,  x,y in R2
%   where r:=x-y,  x,y in R2
%
% target x, unit target normal nx, and source y, are given as complex
% numbers (C \equiv R2). source normal ny unused.
% 
% Mostly useful as reference, & for adaptive integration

% Barnett 6/28/21
G = -StoDLPvelker(mu,y,x,nx)
