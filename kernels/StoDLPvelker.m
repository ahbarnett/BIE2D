function G = StoDLPvelker(mu,x,y,ny)
% STODLPVELKER   kernel func of Stokes double-layer potential, velocity
%
% G = StoDLPvelker(mu,x,y,ny) evaluates
%
%   G(x,y) = (1/pi) (r⊗r) (r.n_y)/ |r|^4 ,  where r:=x-y,  x,y in R2
%
% target x, and source y (with unit normal n_y), are given as complex numbers
% (C \equiv R2). mu is unused.
%
% Mostly useful as reference, & for adaptive integration. Not called in a way
% that needs outer products.

% Barnett 4/26/21
r = x-y;
irr = 1./(conj(r).*r);
rdotny = real(conj(r).*ny);                % r.n_y
rdotnir4 = (1/pi) * rdotny.*irr.*irr;      % w/ overall prefac
G = [real(r).*real(r).*rdotnir4, real(r).*imag(r).*rdotnir4; 
     imag(r).*real(r).*rdotnir4, imag(r).*imag(r).*rdotnir4];
