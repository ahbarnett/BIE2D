function [U L R B T] = doublywalls(U,M)
% DOUBLYWALLS   setup walls with m nodes for doubly-periodic unit cell U
%
% [U L R B T] = doublywalls(U,m) sets up struct U so contains 4 wall segments
%  with m Legendre nodes on each, and normals to right or up. The nodes L.x
%  etc are col vecs but the weights L.w etc row vecs (unlike in setupquad).

% Barnett broken out 6/29/16

[x w] = gauss(M); w=w/2;
L.x = (-U.e1 + U.e2*x)/2; L.nx = (-1i*U.e2)/abs(U.e2) + 0*L.x; L.w=w*abs(U.e2);
R = L; R.x = L.x + U.e1;
B.x = (-U.e2 + U.e1*x)/2; B.nx = (1i*U.e1)/abs(U.e1) + 0*B.x; B.w=w*abs(U.e1);
T = B; T.x = B.x + U.e2;
U.L = L; U.T = T; U.R = R; U.B = B;
