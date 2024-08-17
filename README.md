## PhD Thesis code

This repository contains the MATLAB code for reproducing part of the numerical experiments in my PhD thesis [2].


### Chapter 5 -- A low-memory Lanczos method with rational Krylov compression

The folder `chap5` contains the code of the low-memory Krylov method `RKcompress_fAb` for computing $f(A) \textbf{b}$, which combines an outer Lanczos iteration with an inner rational Krylov subspace used to compress the basis, as well as the code for reproducing the numerical experiments in the chapter. This code was previously used in the preprint [1], and it is also available in the Github repository https://github.com/casulli/ratkrylov-compress-matfun, which also contains a block version.

#### Dependencies

There are no external dependencies for the main function `RKcompress_fAb`. The demo scripts have the following dependencies, required for the comparison against some other low-memory methods:

- `chebfun` https://www.chebfun.org/ (for the AAA algorithm for rational approximation, used in the multishift CG algorithm)
- `funm_quad` http://www.guettel.com/funm_quad/ (for comparison against the funm_quad algorithm)


### Chapter 6 -- Rational Krylov methods for fractional diffusion problems on graphs

The folder `chap6` contains code for solving a fractional diffusion equation on a graph, whose solution is given in the form $f(A) \textbf{b}$, with $f(z) = \exp(-t z^\alpha)$, with $t > 0$ and $\alpha \in $(0,1)$. The code employs desingularization techniques to speed up the convergence of Krylov subspace methods. The folder also includes the scripts for reproducing the figures and tables contained in the chapter. 

#### Dependencies

The code in this folder has the following dependences:

- the `rat_krylov` function from the `RKtoolbox` (http://guettel.com/rktoolbox/) is used to generate the rational Krylov subspaces
- the matrices used in the test problems are available from the SuiteSparse Matrix Collection (https://sparse.tamu.edu/)


### Chapter 9 -- Estimation of gaps in the spectrum

The folder `chap9` contains code for locating gaps in the spectrum of a real symmetric matrix by approximating the traces of spectral projectors associated with sections of its spectrum, using a combination of Hutchinson's stochastic trace estimator and the Lanczos algorithm for computing quadratic forms. The folder contains the main function `gapfinder_main` and the scripts for running all the experiments contained in the chapter.

#### Dependencies

The main algorithm `gapfinder_main` has no external dependencies, but one of the test problems uses a matrix that is available from the SuiteSparse Matrix Collection (https://sparse.tamu.edu/).

### References

[1] Angelo A. Casulli, Igor Simunec, A low-memory Lanczos method with rational Krylov compression for matrix functions, arXiv:2403.04390 (2024).

[2] Igor Simunec, Advances in polynomial and rational Krylov methods for matrix functions with applications, PhD thesis, Scuola Normale Superiore, Italy, 2024.

