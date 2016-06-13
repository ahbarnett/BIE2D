function [u un] = LapDLPeval(t,s,dens)
% LAPDLPEVAL   Evaluate Laplace double-layer potential from curve to targets
%
% [u un] = LapDLPeval(t,s,dens)
%
% Crude native quadr and O(NM) RAM for now
% todo: make O(N+M) & incorporate Gary's scf

% Barnett 6/12/16

if nargout==1
  u = LapDLPmat(t,s) * dens;
else
  [D T] = LapDLPmat(t,s);
  u = D*dens;
  un = T*dens;
end
