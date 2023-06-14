function [ U, d, V ] = svd_sorted( A )
% [ U, d, V ] = svd_sorted( A ) performs an economy-sized SVD of A and
% sorts the singular vectors and singular values in a descending manner.
% 
% IN:
% A   matrix to be decomposed
%
% OUT
% U   matrix of left signular vectors
% d   vector of singular values
% V   matrix of left signular vectors


[U,D,V] = svd(A,'econ');
d = diag(D);
[d, sortIdx] = sort(d, 'descend');
U = U(:, sortIdx);
V = V(:, sortIdx);

end