function [ H_tall, H_fat, C_rel, eps_H ] = lr_approx( H, C_rel )
% [ H_tall, H_fat, C_rel, eps_H ] = lr_approx( H, C_rel ) performs a
% low-rank approximation of H.
%
% IN:
% H         SRP  matrix
% C_rel     array of targeted relative complexities
% 
% OUT: 
% H_tall    cell array, one cell per relative complexity cointaining the
%           corresponding tall component of the LR approximation
% H_fat     cell array, one cell per relative complexity cointaining the
%           corresponding fat component of the LR approximation
% C_rel     actual relative complexities (correspond to target relative
%           compexities up to round-off errors)
% eps_H     matrix approximation error


return_eps_H = (nargout == 4);

% get dimensions
[J, PKmin1] = size(H);

% rank required for targeted C_rel
R = round(C_rel*(2*J*PKmin1)/(2*J + 4*PKmin1));

% actual C_rel
C_rel = (2*J*R + 4*R*PKmin1)/(2*J*PKmin1);

% init
H_tall = cell(length(C_rel),1);
H_fat = cell(length(C_rel),1);
if return_eps_H
    eps_H = nan(length(R), 1);
end

% compute SVD
[U, d, V] = svd_sorted(H);    

% define tall and fat matrix components
for idx = 1:length(R)
    
    H_tall{idx} = U(:, 1:R(idx))*diag(d((1:R(idx))));
    H_fat{idx} = V(:, 1:R(idx))';
    
    % approximation error
    if return_eps_H
        H_approx = H_tall{idx}*H_fat{idx};
        eps_H(idx) = 10*log10(norm(H_approx - H, 'fro').^2/norm(H, 'fro').^2);
    end
    
end

end