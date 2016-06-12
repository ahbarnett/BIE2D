# BIE2D: MATLAB tools for boundary integral equations in 2D

This set of codes solves piecewise constant linear PDEs using boundary integral equations. It includes Laplace, Helmholtz, Stokes, log-singular self-evaluations, and close-evaluation quadratures based on the Cauchy kernel. Also are included periodic BVP solvers. It is somewhat tutorial in nature, and could be used as a template for faster Fortran/C versions.


Code mostly by Alex Barnett, based on work from 2008-2016.

Also includes code by L Nick Trefethen, Gary Marple, Bowei Wu.

## Installation

Download using `git`, `svn`, or as a zip (see green button above).

Open MATLAB in the top level (`BIE2D`) directory, and run `setuppath` to add all needed directories to your path. 

Test by running `fig_cauchycompeval` (for now) which should produce some convergence plots.

Codes have been tested on MATLAB versions from R2012a onwards.

### Directories

`kernels` : Laplace, Stokes, Helmholtz, Cauchy potential evaluators and matrix filling  
`utils`   : general numerical utilities  
`test`    : test codes other than built-in self-tests, figure-generating codes  
`doublyperiodic` : codes for flow (Laplace, Stokes) in doubly-periodic geometries, computation of effective permeability.  
