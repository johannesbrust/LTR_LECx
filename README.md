# LTR_LECx
**Large-Scale Trust-Region Methods For Linear Equality Constrained Optimization**
(extended)

MATLAB codes

This project contains the implementations of the Algorithms from the article
(with possible name changes)

"Large-Scale Optimzation with Linear Equality Constraints", 

J.J. Brust, R.F. Marcia, C.G. Petra, M.A. Saunders 2020/2021.

The folders comprising this project are:

1. examples: Two example scripts to try the algorithms.
1. main: Implementations of Algorithm 1, and Algorithm 2.
    - LTRL2_LEC_V2.m (Alg. 1)
    - LTRSC_LEC_V2.m (Alg. 2)
    - Other solvers for comparisons (LTRL2_LEC,LTRL2_LEC_V1,LTRSC_LEC,LTRSC_LEC_V1)
1. 	auxiliary: General auxiliary routines.
1.	solvers: External Algorithms for comparisons.
1. 	ssget: Interface to the Sparse Suite Matrix collection
        link [SSMat](https://sparse.tamu.edu/).
1.  experiments: 
    - EXPERIMENT_III (Rosenbrock + Suite Sparse A)
    - EXPERIMENT_IV  (CUTEst objective + sprandn A)
1. 	tests: .
1. 	data: Data generated by running the experiments.
1.	figs: Figures generated by the experiments, and used in the article.

## Example
You can run an example from within the /examples folder. On
the MATLAB command prompt:

`> TR1H_QUADRATIC`

```
**********************
Running LTRL2-LEC
it	 obj		 norm(b)	 norm(Proj(g))) 	 norm(dx)	 trrad
0	 2.748e+04	 4.268e-14	 7.296e+03	---		 ---
1	 1.425e+04	 3.404e-14	 3.924e+03	 4.000e+00	 1.544e+01 
2	 1.097e+04	 3.513e-14	 2.149e+03	 1.202e+00	 8.000e+00
3	 8.360e+03	 3.377e-14	 1.512e+03	 2.038e+00	 8.000e+00
4	 6.966e+03	 3.643e-14	 1.323e+03	 1.685e+00	 8.000e+00
5	 5.834e+03	 4.154e-14	 1.607e+03	 2.697e+00	 8.000e+00
6	 5.111e+03	 3.751e-14	 8.800e+02	 1.428e+00	 8.000e+00
7	 4.744e+03	 3.943e-14	 7.296e+02	 6.997e-01	 8.000e+00
8	 4.208e+03	 4.219e-14	 7.947e+02	 2.007e+00	 8.000e+00
9	 4.119e+03	 4.737e-14	 1.360e+03	 1.888e+00	 9.442e-01
10	 3.837e+03	 3.850e-14	 5.000e+02	 4.708e-01	 9.442e-01
11	 3.715e+03	 4.359e-14	 3.623e+02	 4.105e-01	 9.442e-01
12	 3.561e+03	 3.820e-14	 4.508e+02	 9.934e-01	 1.888e+00`
.        .               .               .               .               .
.        .               .               .               .               .
.        .               .               .               .               .
```
