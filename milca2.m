% Mutual Information Least-Dependent Component Analysis (MILCA)
% Reference:
% "Least Dependent Component Analysis Based on Mutual Information"  
% Stogbauer, H.; Kraskov, A.; Astakhov, S. A.; Grassberger, P. 
% Phys. Rev. E, 2004, 70, 066123.
% https://doi.org/10.1103/PhysRevE.70.066123
% https://arxiv.org/abs/physics/0405044

% Author:  Sergey Astakhov (astakhov@gmail.com)
% License: BSD 3 clause

function [Y, W] = milca2(X, threads, P , neighbours, angles, pass)
% Arguments:
% X - an M by N matrix of M mixtures with N observations each 
% threads - the number of CPU threads to use (optional, default = 1), 0 means 
% use all available cores, -1 means all available except 1
% P - the number of components to retain in PCA dimension reduction (optional)
% neighbours - the number of nearest neighbours for mutual information 
% estimator (optional, default 10)
% angles - precision parameter, 2^angles defines the number of angles 
% in MI minimization (optional, default 7)
% pass - smoothing parameter, 2^pass defines the number of lower frequencies 
% to retain in Fourier-approximation of the angular MI landscape (optional, 
% default - half spectrum)
% Output:
% Y - estimated mixture components
% W - de-mixing transformation matrix 
 

    [M,N] = size(X);

    if ~exist('P'), P = M; end
    if ~exist('threads'), threads = 1; end
    if ~exist('neighbours'), neighbours = 10; end
    if ~exist('angles')
        Na = 128;
    else
        Na = 2 ^ angles;
    end
    if ~exist('pass'), pass = int16(Na / 4); end
    
    V = pca_reduction (X, P);
    X = V * X;

    R = minimize_mi (X, Na, pass, neighbours, threads);

    Y = R * X;
    W = R * V;

endfunction


function R = minimize_mi (X, Na, pass, neighbours, threads)
% Minimizes pairwise mutual information of multivariate X over all 2D rotations.
% For parameters Na, pass, neighbours see minimize_mi_2D ()
% Returns cumulative rotation matrix. 

    [P N] = size (X);
    R = eye(P);
    for i = 1:P
      for j = setdiff(1:P, i)
        R2 = eye(P);
        R2([i j], [i j]) = minimize_mi_2D (X([i j], :), Na, pass, neighbours, threads);
        R = R2 * R;
        X = R2 * X;
      endfor
    endfor

endfunction


function R = minimize_mi_2D (X, Na, pass, neighbours, threads)

% Minimizes mutual information of bivariate X over 2D rotation angle.
% The angle is sampled with Na (a power of 2) points first, 
% then the minumum is searched over the smoothed Fourier fit with only
% the lowest pass terms retained (pass has to be less then Na/2).
% neighbours - is the KNN parameter of KSG mutual information estimator.
% Returns rotation matrix.     

  [M N] = size (X);
  delta = 0.5 * pi / (Na - 1);
  Y = zeros(Na * 2, N);
  j = 1;
  for i = 1 : Na
    angle = (i - 1) * delta;
    R = rotation_matrix_2D (angle);
    Y(j:j + 1, :) = R * X;
    j += 2;
  endfor
  
  MI = mi_2D(Y', neighbours, threads);
    
  MI_smooth = low_pass_filter (MI', pass);
  [MI_min min_indx] = min (MI_smooth);
  angle = (min_indx - 1) * delta;
  
  R = rotation_matrix_2D (angle);
 
endfunction


function R = rotation_matrix_2D (angle)
  
  R = [cos(angle) sin(angle); -sin(angle) cos(angle)];  

endfunction


function y = low_pass_filter (x, pass)
% Smoothens a univariate signal by a low-pass Fourier fit.
% x - is assumed to be a (1:N) vector 
% pass - the number of lower FFT bins to retain, 
% higher frequencies get suppressed
  
  N = size(x, 2);
  s = fft(x, N);

  window = zeros(1, N);
  window(1:pass + 1) = 1; 
  window(end - pass + 1:end) = 1;     
  
  y = ifft(s.*window, N);
  
endfunction


