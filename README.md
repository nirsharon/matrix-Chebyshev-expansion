# matrix-Chebyshev-expansion
software for the paper: "Matrix Chebyshev expansion and its application for eigenspaces recovery" 
  
Authors:               Nir Sharon
                            Yoel Shkolnisky
                            
Version: 1.0.

Release date: 31 August 2016.

Matlab code.

-------------------------------------------------------------------
 Usage
-------------------------------------------------------------------
There are 9 scripts that generate figures and table of the paper. Run each separately, change “saveit” variable to 0 to ensure no new files are saved.

Other scripts support the main scripts by generating the Chebyshev coefficients, evaluating the matrix expansion or evaluating the vector-matrix version of the latter.

-------------------------------------------------------------------
 3rd party code
-------------------------------------------------------------------
This package uses and is distributed the following software produced by 3rd parties:
1. Chebyshev coefficients evaluation from the Numerical Recipes book.

-------------------------------------------------------------------
main scripts descriptions 
-------------------------------------------------------------------
non_smooth_f1                 — Figure 1, part1: the first non smooth function f_1
non_smooth_f2                 — Figure 1, part2: the first non smooth function f_2
non_smooth_f2                 — Figure 2, both parts: numerics on function f_3
JordanDecay2          	      — Figure 3, part 1: Jordan blocks of different sizes
JordanDecay_differentEV — Figure 3, part2: Jordan blocks of different eigenvalues.
MultiplyBlocks                   — Figure 4
plot_FilterFunc         		   —  Figure 5, the plot of the filter function
EIGvsChebyshev     	      —  Figure 6, the comparison of Chebyshev and EIG
HugeMatrixExample 		   —  generates Table 1
