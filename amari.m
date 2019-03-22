% Computes the Amari performance index (used for assessing mixture
% decomposition quality on the ground truth and estimated mixing matrices). 
% It is a "distance" measure between 2 matrices which is insentive 
% to permutation and scaling.
% Empirically, for a successful decomposition the Amari index
% falls well below 0.1

% Reference:
% S. Amari, A. Cichocki and H.H. Yang, A new learning
% algorithm for blind source separation, in D.S. Touretzky
% et al. eds., Advances in Neural Information Processing
% 8 (Proc. NIPS’95), pp. 757-763 (MIT Press, Cambridge,MA, 1996).

% Author:  Sergey Astakhov (astakhov@gmail.com)
% License: BSD 3 clause

function out = amari(C, A)

    [b, a] =  size(C);

    d1 = pinv(A) * C;
    d1 = sum(ntu(abs(d1))) - 1;

    d2 = pinv(C) * A;
    d2 = sum(ntu(abs(d2)))-1;

    out = (sum(d1) + sum(d2)) / (2 * a * (a - 1));

endfunction

function CN = ntu(C)

    [m n] = size(C);
    for t=1:n
     CN(:, t) = C(:, t)./max(abs(C(:, t)));
    end

endfunction