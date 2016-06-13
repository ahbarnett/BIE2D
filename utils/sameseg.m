function a = sameseg(t,s)
% SAMEREG   Return true if two segments have numerically the same nodes

% Barnett 6/12/16
a = (numel(s.x)==numel(t.x)) && max(abs(s.x-t.x))<1e-14;
