function setuppath
% SETUPPATH   puts all M-files from BIE2D package in path as absolute paths
%
% Barnett 6/12/16
mfilepath=fileparts(mfilename('fullpath'));
addpath([mfilepath, '/kernels']);
addpath([mfilepath, '/utils']);
addpath([mfilepath, '/test']);
addpath([mfilepath, '/solvers']);
addpath([mfilepath, '/doublyperiodic']);
