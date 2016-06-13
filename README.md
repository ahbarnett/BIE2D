# BIE2D: MATLAB tools for boundary integral equations on curves in 2D

This set of codes solves piecewise constant linear PDEs using boundary integral equations on curves.  It includes Laplace, Helmholtz, Stokes, log-singular self-evaluations, and close-evaluation quadratures based on the Cauchy kernel.  Also are included periodic BVP solvers.  It is designed to be reasonably efficient yet tutorial in nature (easy to read and well-documented), and could be used as a template for faster Fortran/C versions.


Code mostly by Alex Barnett, based on work from 2008-2016 
including MPSpack, LSC2D and other projects.

Also includes the following contributions:

  Gary Marple - matrix versions of global close evaluation quadratures  
  Bowei Wu - Stokes velocity extension from Laplace  
  L Nick Trefethen - Gaussian quadrature  


## Installation

Download using `git`, `svn`, or as a zip (see green button above).

Open MATLAB in the top level (`BIE2D`) directory, and run `setuppath` to add all needed directories to your path. 

Test by running `testall` which should produce lots of outputs and convergence plots without crashing.

Codes have been tested on MATLAB versions from R2012a onwards.


## Directories

`kernels` : Laplace, Stokes, Helmholtz, Cauchy potential evaluators and matrix filling  
`utils`   : general numerical utilities  
`test`    : test codes other than built-in self-tests, figure-generating codes
`solvers` : 2D BVP solver example codes, also serve to test kernels    
`doublyperiodic` : codes for flow (Laplace, Stokes) in doubly-periodic geometries, computation of effective permeability.  


## Improvements needed

* Helmholtz  
* bring in QBX?  
