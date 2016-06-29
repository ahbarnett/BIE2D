function s = wobblycurve(a,w,N)
% WOBBLYCURVE   Set up a wobbly smooth closed curve ("starfish")
%
% s = wobblycurve(a,w) where a is amplitude of wobble (eg 0.3) and w is
%  the frequency (eg 5), returns a segment struct of the form as in
%  setupquad, but with s.inside a handle to Boolean test whether a point is
%  inside the curve.
%
% Without arguments, does self-test of analytic 1st and 2nd derivatives

% Barnett repackaged 6/12/16, generic angle offset 6/29/16
if nargin==0, test_wobblycurve; return; end

th = 0.2;    % generic rotation
% since analytic derivs not too messy, use them...
R = @(t) 1 + a*cos(w*(t-th));
Rp = @(t) -w*a*sin(w*(t-th));
Rpp = @(t) -w*w*a*cos(w*(t-th));
s.Z = @(t) R(t).*exp(1i*t);
s.Zp = @(t) (Rp(t) + 1i*R(t)).*exp(1i*t);
s.Zpp = @(t) (Rpp(t) + 2i*Rp(t) - R(t)).*exp(1i*t);
s = setupquad(s,N);
s.inside = @(z) abs(z)<R(angle(z));

%%%%%%%
function test_wobblycurve
s = wobblycurve(0.3,5,150);
max(abs(s.xp - perispecdiff(s.x)))     % make sure analytic close to numerical
max(abs(s.xpp - perispecdiff(s.xp)))
