// Mutual information estimators 

// Author:  Sergey Astakhov (astakhov@gmail.com)
// License: BSD 3 clause

#include <stdlib.h>
#include <float.h>

#include "miutils.h"

double compute_mi_2D_ksg(double *data, int N, int K)
// Reference:
// "Estimating Mutual Information" 
// A. Kraskov, H. Stogbauer and P. Grassberger, 
// Phys. Rev. E 69, 066138 (2004). 
// https://arxiv.org/abs/cond-mat/0305641
// https://doi.org/10.1103/PhysRevE.69.066138

{
  double **x, *psi, mi, scalxx[2];
  
  x = (double**) malloc(2 * sizeof(double*));
  x[0]=(double*) malloc(N * sizeof(double));
  x[1]=(double*) malloc(N * sizeof(double));
  
  for (int i = 0; i < N; ++i) 
  {
       x[0][i] = data[i];
       x[1][i] = data[i + N];
  } 
        
  scalxx[0] = scalxx[1] = - (N - 5)/DBL_MAX;
  psi=(double*) calloc(N + 1, sizeof(double)); 
  psi[1] = -(double).57721566490153;
  for (int i = 1; i < N; i++) psi[i + 1] = psi[i] + 1 / (double)i;

  mi2r(x, N, K, psi, &scalxx[0], &mi);
  
  free(x[0]); free(x[1]); free(x); free(psi);
    
  return mi;
}
