function [ varargout ] = slci_approx( H, Delta_rel, N_p, C_rel, constraint )
% [ varargout ] = slci_approx( H, Delta_rel, N_p, C_rel, constraint )
% performs a sampling + low-complexity interpolation approximation of H,
% where low complexity is achieved by either low-rank or sparse
% interpolation.
%
% IN:
% H             SRP  matrix
% Delta_rel     TDOAs relative to sampling period - candidate locations x 
%               mic. pairs
% N_p           number of required samples per mic. pair
% C_rel         array of targeted relative complexities
% constraint    'LR' or 'SP', indicating low-rank or sparse approximation
%               of the interpolation matrix
%
% OUT: 
% Lambda_tall   tall component of interpolation matrix - candidate
%               locations x rank (available if constraint = 'LR')
% Lambda_fat    tall component of interpolation matrix - rank x total
%               number of samples (available if constraint = 'LR')       
% Lambda_sp     sparse interpolation matrix - candidate locations x total
%               number of samples (available if constraint = 'SP') 
% S             sampling matrix
% C_rel         actual relative complexities (correspond to target relative
%               compexities up to round-off errors)
% eps_H         matrix approximation error




% get dimensions
[J, PKmin1] = size(H);
P = size(Delta_rel, 2);
K = PKmin1/P + 1;

% perform SI approximation
[ Lambda, S ] = si_approx( H, Delta_rel, N_p );

 % total number of required samples
N = size(S,1);             

% complexity of sampling (assuming use of iFFT)
C_samp = 8*P*K*log2(2*K);

% init Lambda output and get R or Q
switch constraint
    
    case 'LR'
        
        return_eps_H = (nargout >= 5);
        
        Lambda_tall = cell(length(C_rel),1);
        Lambda_fat = cell(length(C_rel),1);
        
        % get R
        R = round((C_rel*(2*J*PKmin1) - C_samp)/(J + N));
        
        % actual C_rel
        C_rel = (R*(J + N) + C_samp)/(2*J*PKmin1);
        
    case 'SP'
        
        return_eps_H = (nargout >= 4);
        
        Lambda_sp = cell(length(C_rel),1);
        
        % get Q
        Q = round(C_rel*(2*J*PKmin1) - C_samp);
        
        % actual C_rel
        C_rel = (Q + C_samp)/(2*J*PKmin1);
        
end

% not all R or Q are reasonable
startIdx = find(C_rel > C_samp/(2*J*PKmin1), 1, 'first');            % lower bounded by complexity of samp
endIdx   = find(C_rel < (J*N + C_samp)/(2*J*PKmin1), 1, 'last');     % upper bounded by complexity of samp + inter

if isempty(startIdx)
    error('all relative complexities below lower bound.')
end
if isempty(endIdx)
    error('all relative complexities above upper bound.')
end

% use only reasonable R and Q
switch constraint
    case 'LR'
        R(1:startIdx-1) = nan;
        R(endIdx+1:end) = nan;
    case 'SP'
        Q(1:startIdx-1) = nan;
        Q(endIdx+1:end) = nan;
end

% init error
if return_eps_H
    eps_H = nan(length(C_rel),1);
end

% solve min problem
switch constraint
    case 'LR'
        
        % perform SVD
        [U,d,V] = svd_sorted(Lambda);
        
        % build Lambda_tall and Lambda_fat
        for idx = startIdx:endIdx
            
            Lambda_tall{idx} = U(:, 1:R(idx))*diag(d((1:R(idx))));
            Lambda_fat{idx} = V(:, 1:R(idx))';
            
            % approximation error
            if return_eps_H
                Lambda_lc = Lambda_tall{idx}*Lambda_fat{idx};
                H_approx =  Lambda_lc*S;
                eps_H(idx) = 10*log10(norm(H_approx - H, 'fro').^2/norm(H, 'fro').^2);
            end
            
        end
        
    case 'SP'

        % sort elements
        Lambda_vec = Lambda(:);
        [~, sortIdx] = sort(abs(Lambda_vec), 'descend');

        
        % build Lambda_sp
        for idx = startIdx:endIdx
            Lambda_vec_sp = Lambda_vec;
            Lambda_vec_sp(sortIdx(Q(idx)+1:end)) = 0;
            Lambda_sp{idx} = sparse(reshape(Lambda_vec_sp, size(Lambda)));
            % approximation error
            if return_eps_H
                Lambda_lc = Lambda_sp{idx};
                H_approx = full(Lambda_lc)*S;       % in MATLAB, Lambda_lc*S with Lambda_lc being of type sparse double is inefficient if Lambda_lc is not very sparse. As we also test higher complexities, we use a full matrix here instead. 
                eps_H(idx) = 10*log10(norm(H_approx - H, 'fro').^2/norm(H, 'fro').^2);
            end
        end
        
end

% wite output variables
switch constraint
    case 'LR'
        varargout{1} = Lambda_tall;
        varargout{2} = Lambda_fat;
        varargout{3} = S;
        varargout{4} = C_rel;
        if return_eps_H
            varargout{5} = eps_H;
        end
    case 'SP'
        varargout{1} = Lambda_sp;
        varargout{2} = S;
        varargout{3} = C_rel;
        if return_eps_H
            varargout{4} = eps_H;
        end
end

end