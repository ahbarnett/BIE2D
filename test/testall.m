% complete set of tests for BIE2D. All should produce small numbers &
% occasional figures and tables.

% utils
perispecdiff
perispecinterp
setupquad
showsegment
wobblycurve

% kernels: basics (native and self rules)
testGRFLapHelm
testStokernels

% kernels: close eval
Cau_closeglobal
fig_Cau_closeglobal
LapSLP_closeglobal
LapDLP_closeglobal
StoSLP_closeglobal
StoDLP_closeglobal

% solvers: native quadr and close eval
LapintDirBVP
StointDirBVP

