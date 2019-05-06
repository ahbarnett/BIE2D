# BIE2D: MATLAB tools for boundary integral equations on curves in 2D

This set of codes solves boundary value problems for piecewise constant coefficient linear PDEs using potential theory, ie boundary integral equations (BIE) on curves. Quadratures that are very high-order or spectral are used, allowing accuracies approaching machine precision with small numbers of unknowns. The idea is to provide a simple and uniform interface to the evaluation of layer potentials and filling of Nystrom matrices for Laplace, Helmholtz, and Stokes kernels, including modern quadratures for singular self-evaluation, and close-evaluation quadratures (eg based on the Cauchy kernel).  Simple BVP solvers are included, for various geometries including singly and doubly periodic. This provides a sandbox to deliver implementations of various schemes that are under active research by myself and collaborators, as well as by other experts such as R. Kress, J. Helsing. The code is designed to be reasonably efficient, yet tutorial in nature (easy to read and well-documented), and could be used as a template for faster Fortran/C versions.  I plan to include fast algorithms, panel quadratures (including close-evaluations), and corner quadratures.

Main author: Alex Barnett

Version: 20190506

Based on work from 2008-2016 including subsuming the integral-equation parts of [MPSpack](https://github.com/ahbarnett/mpspack), all of [LSC2D](http://math.dartmouth.edu/~ahb/software/lsc2d.tgz), and my [BIE tutorial](https://math.dartmouth.edu/~fastdirect/notes/quadrtut.zip).

Also includes the following contributions and influences:

  Bowei Wu - Stokes velocity extension from Laplace  
  Gary Marple - matrix versions of global close evaluation quadratures  
  Nick Trefethen - Gaussian quadrature
  Jun Wang - 2nd derivs of Cauchy, Laplace SLP, and traction of Stokes SLP  

As of 2018, David Stein made a python implementation of most of BIE2D, including
new fast versions of close evaluations, in his package [pyBIE2D](https://github.com/dbstein/pybie2d).


### Installation and testing

Download using `git`, `svn`, or as a zip (see green button above).

Open MATLAB in the top level (`BIE2D`) directory, and run `bie2dsetup` to add all needed directories to your path. 

Test by running `testall` which should produce lots of error outputs close to machine precision, figures, etc, and yet not crash.

Many functions (eg close routines in kernels) have built-in self-tests run by calling without any arguments.

Codes have not been tested on MATLAB versions prior to R2012a.


### Directories

`kernels` : Laplace, Stokes, Cauchy potential evaluation and matrix filling, including close-evaluation  
`utils`   : general numerical and plot utilities  
`test`    : test codes (other than built-in self-tests), figure-generating codes  
`solvers` : 2D BVP solver example codes, also serve to test kernels  
`solvers/closetouchingexpts` : close-to-touching curve experiments  
`singlyperiodic` : codes for Laplace, Stokes flow in periodic pipes, possibly with vesicles (to do)  
`doublyperiodic` : codes for flow (Laplace, Stokes) in doubly-periodic geometries, computation of effective permeability (in progress)  

### Notes and design decisions

1. Every kernel can be accessed as a dense matrix mapping density to potential (if no density function is given), or as the potential evaluation given a single density (or stack of density columns).  These are bundled into the same calling interface.  Usually one of these simply calls the other (eg the matrix filling for close-evaluation is currently done by sending in columns of the identity as densities, which is wasteful), but I plan to provide more efficient matrix fillers based on BLAS3 later.

1. For all coordinates in $\mathbb{R}^2$ we use complex numbers in the form $x+iy$, since this is very convenient for geometry. For Cauchy integral interpolation we obviously also use complex numbers.

1. Stokes involves vector-valued densities, velocities, etc. For clarity of coding and visualization we have decided to order the node index fast and the vector index slow, ie to stack all the x-components, followed by all the y-components. Matrices thus have a block structure of the form $[A_{11}, A_{12}; A_{21}, A_{22}]$. The other choice of alternating x and y components would be better for RAM locality, and to feed into direct solvers, but we err instead on the side of simplicity/reability and leave this for a Fortran/C implementation. Likewise, we have avoided the use of complex numbers to represent x and y components for Stokes (apart from in internal routines such as Laplace extension).

1. So far everything is based on the global periodic trapezoid rule. Thinking ahead, will panels be cell arrays of segments, or one large segment struct with breakpoints listed?

### Action items

* decide where to build in switches for close eval - per target? (scf)
* multiple closed curves (islands) helpers, use for dpls figs?
* Green's representation theorem kernel tests for Lap, then Sto (ugh), so don't rely on BIE density soln for basic kernel tests
* more BVP solver demos (eg bring over testStokesSDevalclose.m w/ all 4 BVPs)
* FMM MEX interfaces
* [long-term] Convert whole thing to C/Fortran libraries
* [long term] basic fast direct solver examples - jj index list fields in s,t?
* kd-tree for close-evaluation lists (Marple). Need non-Toolbox kd-tree.
* [low-priority] Alpert and other options for log-singular kernels
* Helmholtz bring in from MPSpack
* panels (Helsing-style), corners with panels, bring in from various tests
* [long term] MEX interface to Rachh QBX/FMM ?
* [low-priority] bring in singlyperiodic, get Adrianna codes

### Done

* Cau_closeglobal simpler uses "interpolate the derivative", exterior S-W form
* cleaner kernel interface without mat or eval suffices
* many repmats/ones -> bsxfun for speed in kernels
* derivs for 'e' LapSLP_closeglobal corrected for nonzero totchg.
* Stokes mat fills bsxfun, StoSLP, StoSLP_closeglobal debug via vel BVP for now
* Stokes SLP pressure closeglobal eval, upsampling
* StoDLP_cg self-test & pressure close, StointvelBVP all pres plots added
* finish bring over doublyperiodic, pres figs
* srcsum2 for speed (targ sum not src, better for close eval)
* initial close-touching expts
* much accelerated the Cauchy close global matrix fill, using BLAS3 for all the O(N^3) parts, hence accelerating Laplace & Stokes (which call Cauchy in non-sparse way due to CSLP matrix being dense)
