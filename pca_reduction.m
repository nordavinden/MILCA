% PCA with dimension reduction

% Author:  Sergey Astakhov (astakhov@gmail.com)
% License: BSD 3 clause

function V = pca_reduction (X, P)
%Performs principal component analysis on X with dimension reduction
%into P components with largest covariance eigenvalues
  
  [M, N] = size(X); 
  [E_, D_] = eig(cov(X'));
  V = inv(sqrtm(D_((M - P + 1):M, (M - P + 1):M))) * E_(:, (M - P + 1):M)';
  
endfunction
