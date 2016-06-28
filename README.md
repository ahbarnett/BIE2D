# BIE2D: MATLAB tools for boundary integral equations on curves in 2D

This set of codes solves boundary value problems for piecewise constant coefficient linear PDEs using potential theory, ie boundary integral equations (BIE) on curves. Quadratures that are very high-order or spectral are used, allowing accuracies approaching machine precision with small numbers of unknowns. The idea is to provide a simple and uniform interface to the evaluation of layer potentials and filling of Nystrom matrices for Laplace, Helmholtz, and Stokes kernels, including modern quadratures for singular self-evaluation, and close-evaluation quadratures (eg based on the Cauchy kernel).  Simple BVP solvers are included, for various geometries including singly and doubly periodic. This provides a sandbox to deliver implementations of various schemes that are under active research by myself and collaborators, as well as by other experts such as R. Kress, J. Helsing. The code is designed to be reasonably efficient, yet tutorial in nature (easy to read and well-documented), and could be used as a template for faster Fortran/C versions.  I plan to include fast algorithsm, panel quadratures (including close-evaluations), and corner quadratures.

Main author: Alex Barnett, based on work from 2008-2016 including subsuming the integral-equation parts of [MPSpack](https://github.com/ahbarnett/mpspack), all of [LSC2D](http://math.dartmouth.edu/~ahb/software/lsc2d.tgz), and my [BIE tutorial](https://math.dartmouth.edu/~fastdirect/notes/quadrtut.zip).

Version: 20160627.

Also includes the following contributions:

  Gary Marple - matrix versions of global close evaluation quadratures  
  Bowei Wu - Stokes velocity extension from Laplace  
  L Nick Trefethen - Gaussian quadrature  

### Installation

Download using `git`, `svn`, or as a zip (see green button above).

Open MATLAB in the top level (`BIE2D`) directory, and run `bie2dsetup` to add all needed directories to your path. 

Test by running `testall` which should produce lots of error outputs close to machine precision, figures, etc, and yet not crash.

Codes have been tested on MATLAB versions from R2012a onwards.


### Directories

`kernels` : Laplace, Stokes, Helmholtz, Cauchy potential evaluators and matrix filling  
`utils`   : general numerical utilities  
`test`    : test codes other than built-in self-tests, figure-generating codes  
`solvers` : 2D BVP solver example codes, also serve to test kernels  
`singlyperiodic` : codes for Laplace, Stokes flow in periodic pipes, possibly with veiscles.  
`doublyperiodic` : codes for flow (Laplace, Stokes) in doubly-periodic geometries, computation of effective permeability.  

### Notes and design decisions

1. Every kernel can be accessed as a dense matrix (if no density function is given), or as the evaluation given a density.  These are bundled into the same calling interface.  Usually one of these simply calls the other, but this allows dropping in more efficient versions for either.

1. For all coordinates in $\mathbb{R}^2$ we use complex numbers in the form $x+iy$, since this is very convenient for geometry. For Cauchy integral interpolation we obviously also use complex numbers.

1. Stokes involves vector-valued densities, velocities, etc. For clarity of coding and visualization we have decided to order the node index fast and the vector index slow, ie to stack all the x-components, followed by all the y-components. Matrices thus have a block structure of the form $[A_{11}, A_{12}; A_{21}, A_{22}]$. The other choice of alternating x and y components would be better for RAM locality, and to feed into direct solvers, but we err instead on the side of simplicity/reability and leave this for a Fortran/C implementation. Likewise, we have avoided the use of complex numbers to represent x and y components for Stokes.

1. Will panels be cell arrays of segments, or one large segment struct with breakpoints listed?

### Improvements needed

* exploit incoming 0s for efficient Cau_closeglobal matrix filling, Lap, Sto too, needing special matrix-filling all the way up to avoid O(MN^2)
* Stokes mat fill bsxfun
* Green's representation theorem tests
* more BVP solver demos
* FMM MEX interfaces
* kd-tree for close-evaluation lists (Marple)
* panels
* Alpert and other options for log-singular kernels
* Helmholtz bring in from MPSpack
* corners with panels, bring in from various tests
* bring in QBX ?

### Done

* Cau_closeglobal simpler uses "interpolate the derivative", exterior S-W form
* cleaner kernel interface without mat or eval suffices
* some repmats -> ones for speed
* derivs for 'e' LapSLP_closeglobal corrected for nonzero totchg.
