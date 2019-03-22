// Propagate computation of pairwise mutual information 
// over a data matrix in a serial or parallel manner

// Author:  Sergey Astakhov (astakhov@gmail.com)
// License: BSD 3 clause

#include <thread>
#include "compute_mi.hpp"

void compute_mi_serial (double *x, double *mi, int M, int N, int K)
{
    for (int j = 0; j < N; j += 2) 
        mi[j / 2] = compute_mi_2D_ksg(&x[j * M], M, K);
}

void compute_mi_parallel (double *x, double *mi, int M, int N, int K, int n_threads)
{
    std::thread t[n_threads];
    int N_res = 2 * ((N / 2) % n_threads);
    int N_part =  2 * (((N - N_res) / 2) / n_threads);
    
    int i, k;
    for (k = 0; k < n_threads; ++k)
    {
        i = N_part * k;
        t[k] = std::thread(compute_mi_serial, 
                           &x[i * M], &mi[i / 2], M, N_part, K);
    }
    i = N_part * k;
    if (i < N) compute_mi_serial(&x[i * M], &mi[i / 2], M, N_res, K);
    
    for (k = 0; k < n_threads; ++k) t[k].join();    
}