function [u un] = LapSLPeval(t,s,dens)
% LAPSLPEVAL   Evaluate Laplace single-layer potential from curve to targets
%
% [u un] = LapSLPeval(t,s,dens)
%
% Crude native quad and O(NM) RAM for now
% todo: make O(N+M) & incorporate Gary's scf

% Barnett 6/12/16

if nargout==1
  u = LapSLPmat(t,s) * dens;
else
  [S DT] = LapSLPmat(t,s);
  u = S*dens;
  un = DT*dens;
end
