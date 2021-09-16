function [pa tpan s] = quadr_uniform_panels(s,Np,p)
% QUADR_UNIFORM_PANELS   setup uniform-in-parameter panel quadrature
%
% [pa tpan s] = quadr_uniform_panels(s,Np,p) where s is a segment struct
% defining a global curve by s.Z and s.Zp over 2pi-periodic parameter,
% with Np panels each of order p,
% returns a struct array pa of panels (each a usual BIE2D segment struct),
% tpan a list of parameter endpoints (length Np+1), and s a modified single
% segment struct with union of all quadr point arrays.
%
% Just a crude helper for now, tested in testGRFLap_panels

% Barnett 9/15/21
[stdpan.t, stdpan.w] = gauss(p);
tpan = linspace(0,2*pi,Np+1);   % panel parameter endpoints
s.x = []; s.xp = []; s.w = []; s.nx = []; s.wxp = [];   % override s nodes
for i=1:Np                    % build each panel, oh, also and quadr info in s
  hl = (tpan(i+1)-tpan(i))/2; % param half-size
  t = tpan(i) + hl*(1+stdpan.t);  % global param vals for this pan
  pa{i}.x = s.Z(t); s.x = [s.x; pa{i}.x];
  pa{i}.xp = s.Zp(t); s.xp = [s.xp; pa{i}.xp];
  pa{i}.nx = -1i*pa{i}.xp./abs(pa{i}.xp); s.nx = [s.nx; pa{i}.nx];
  pa{i}.sp = abs(pa{i}.xp);   % speeds (wrt global param)
  pa{i}.w = hl*stdpan.w.'.*pa{i}.sp; s.w = [s.w; pa{i}.w];  % col vec, speed wei
  pa{i}.wxp = hl*stdpan.w.'.*pa{i}.xp; s.wxp = [s.wxp; pa{i}.wxp]; % Helsing wzp
end
