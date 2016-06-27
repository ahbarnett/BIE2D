% complete set of tests for BIE2D. All should produce small numbers &
% occasional figures and tables.
% Barnett 6/12/16

% utils
perispecdiff
setupquad
showsegment
wobblycurve

% kernels: native quadr
LapintDirBVP
testStokernels
StointDirBVP

% kernels: close eval
Cau_closeglobal
fig_Cau_closeglobal
LapDLP_closeglobal
