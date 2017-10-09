function s = reparam_bunched(s,be)
% REPARAM_BUNCHED  Reparameterize a segment slowing down at 0,pi/2,pi,3pi/2.
%
% s = reparam_bunched(s,be) takes a segment struct s and returns another,
%  where be is the beta parameter giving the angular range devoted to the
%  central. The bunching factor is of order exp(be), or bunching region
%  of order exp(-be).
%
% Note: s.cur may be inaccurate. The user should consider replacing with
%  analytic values
%
% Without arguments, does self-test
%
% Barnett 9/28/17

if nargin==0, test_reparam_bunched; return; end

% set up reparam func on grid in [0,2pi)
t = s.t; N = numel(t); %(1:N)/N * 2*pi;
h = 2*pi/N;
yp = cosh(be*sin(2*t));
I = sum(yp)*h;
al = 2*pi/I;
yp = yp*al;   % make yp integrate to 2*pi
y = perispecint(yp);
y = y-y(1);   % start at t=0

s.x = s.Z(y);     % new nodes and their derivs, exactly
s.xp = yp.*s.Zp(y);
% (here could replace s.xpp analytically)

scopy = s;
s = rmfield(s,{'Z','Zp','Zpp'});  % kill otherwise get used
s = setupquad(s);  % regenerate all from from s.x and s.xp
%s.t = y;           % in case user needs...
%s.Z = scopy.Z;
%s.Zp = scopy.Zp;
%s.Zpp = scopy.Zpp;

%%%%%%%%
function test_reparam_bunched
N=200;
figure;      % animate...
for be=1:0.5:8
  s = wobblycurve(0.5,0,1,N);
  s = reparam_bunched(s,be);
  clf; showsegment(s); title(sprintf('\\beta = %.3f\n',be)); drawnow
end
%s.cur % round-off deviates from 2 by up to 1e-7
