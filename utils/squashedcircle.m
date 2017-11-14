function s = squashedcircle(R,al,b,N)
% SQUASHEDCIRCLE  Set up a rotated squashed circle smooth closed curve.
%
% s = squashedcircle(R,al,b,N) returns a segment structure s.
% Inputs:
% R = max radius
% al = rotation angle
% b in (0,1] is squashing param (0 -> touching)
% N is # bdry nodes
%
% Barnett 11/14/17
if nargin==0, test_squashedcircle; return; end

R = R*exp(1i*al);  % orientation
s.Z = @(t) R*(exp(1i*t) - 1i*(1-b)*sin(t).^3);
s.Zp = @(t) R*(1i*exp(1i*t) - 3i*(1-b)*cos(t).*sin(t).^2);
s.Zpp = @(t) R*(-exp(1i*t) - 3i*(1-b)*(-sin(t).^3 + 2*cos(t).^2.*sin(t)));
s = setupquad(s,N);
s.inside = @(z) abs(imag(z/R)) < sqrt(1-real(z/R).^2) - (1-b)*sqrt(1-real(z/R).^2).^3;

%%%%%%%
function test_squashedcircle
s = squashedcircle(1.5,1.2,0.1,150);
max(abs(s.xp - perispecdiff(s.x)))     % make sure analytic close to numerical
max(abs(s.xpp - perispecdiff(s.xp)))
g = -2:0.01:2; [xx yy]=meshgrid(g); zz=xx+1i*yy;
figure; imagesc(g,g,s.inside(zz)); hold on; plot(s.x,'w.'); axis xy equal
