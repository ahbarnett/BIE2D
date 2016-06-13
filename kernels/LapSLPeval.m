function [u un] = LapSLPeval(t,s,dens)
% LAPSLPEVAL   Evaluate Laplace single-layer potential from curve to targets
%
% [u un] = LapSLPeval(t,s,dens)
%
% t.x is column vector of targets.
% dens is a column vector, or stack of n column vectors.
%
% Crude native quad and O(NM) RAM for now
%
% Tested by: LAPINTDIRBVP

% todo: make O(N+M) & incorporate Gary's scf
% * couple to FMM for each col. * doc formula

% Barnett 6/12/16

if nargout==1
  u = LapSLPmat(t,s) * dens;
else
  [S DT] = LapSLPmat(t,s);
  u = S*dens;
  un = DT*dens;
end
