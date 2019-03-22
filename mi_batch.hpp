// Propagate computation of pairwise mutual information 
// over a data matrix in a serial or parallel manner

// Author:  Sergey Astakhov (astakhov@gmail.com)
// License: BSD 3 clause

void compute_mi_serial (double *x, double *mi, int M, int N, int K);

void compute_mi_parallel (double *x, double *mi, int M, int N, int K, int n_threads);
