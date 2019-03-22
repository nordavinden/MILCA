// A MATLAB/Octave mex interface to batch parallel KNN estimators 
// of mutual information for bivariate data

// Author:  Sergey Astakhov (astakhov@gmail.com)
// License: BSD 3 clause


#include <thread>
#include "mex.h"
#include "mi_batch.hpp"

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Takes 3 arguments
// prhs[0] - data matrix of bivariate signals (N paires of columns) with M rows
// prhs[1] - the number of nearest neighbours for KNN MI estimators (optional, default = 10)
// prhs[2] - the number of CPU threads to use (optional, default = 1), 0 means use all 
// available cores, -1 means all available except 1
// Returns a vector of N mutual information values
{
  
  int M, N, n_neighbours, n_threads;
  double *data, *mi;
  
  if (nrhs == 0) 
    mexErrMsgTxt("At least a data matrix is required in arguments.");
  
  data = (double *)mxGetPr(prhs[0]);
  
  if (nrhs < 3)
  {
      n_threads = 1;
      if (nrhs < 2)
         n_neighbours = 10;
      else 
         n_neighbours = mxGetScalar(prhs[1]);
  }
  else 
  {
      n_threads = mxGetScalar(prhs[2]);
      n_neighbours = mxGetScalar(prhs[1]);
  }
    
  if (n_threads < -1) 
    mexErrMsgTxt("The number of threads needs to be -1, 0 or a positive number.");
  if (n_neighbours <= 0) 
    mexErrMsgTxt("The number of nearest neighbours needs to be a positive number.");

  int threadpool = std::thread::hardware_concurrency();
  if (n_threads == 0) n_threads = threadpool;
  if (n_threads == -1) n_threads = threadpool - 1;
    
  M = mxGetM(prhs[0]); // rows 
  N = mxGetN(prhs[0]); // columns  

  if (N % 2 != 0) 
    mexErrMsgTxt("The input matrix is supposed to have an even number of columns.");

  plhs[0] = mxCreateDoubleMatrix(N / 2, 1, mxREAL);
  mi = mxGetPr(plhs[0]);
      
  if (n_threads == 1 || N == 2) 
      compute_mi_serial (data, mi, M, N, n_neighbours);
  else 
      compute_mi_parallel (data, mi, M, N, n_neighbours, n_threads);
} 
