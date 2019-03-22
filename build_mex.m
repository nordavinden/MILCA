% Build mex routines with their dependencies
% Note that miutils.C is currently not in this repo

% Author:  Sergey Astakhov (astakhov@gmail.com)
% License: BSD 3 clause

clear -f
mex mi_2D.cpp mi_batch.cpp compute_mi.cpp miutils.C -O2