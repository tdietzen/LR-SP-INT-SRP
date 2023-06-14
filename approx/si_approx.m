function [ Lambda, S, C_rel, eps_H ] = si_approx( H, Delta_rel, N_p )
% [ Lambda, S, C_rel, eps_H ] = si_approx( H, Delta_rel, N_p ) performs a
% sampling+interpolation approximation of H.
%
% IN:
% H           SRP  matrix
% Delta_rel   TDOAs relative to sampling period - candidate locations x 
%             mic. pairs
% N_p         number of required samples per mic. pair
% 
% OUT: 
% Lambda      interpolation matrix - candidate locations x total number of 
%             samples
% S           sampling matrix
% C_rel       actual relative complexities (correspond to target relative
%             compexities up to round-off errors)
% eps_H       matrix approximation error

return_eps_H = (nargout == 4);

% get dimensions
J = size(H, 1);
P = size(Delta_rel, 2);
PKmin1 = size(H, 2);
Kmin1 = PKmin1/P;
K = Kmin1+1;

% frequency bin index
k = 1:Kmin1;

% init
S = [];
Lambda = [];

for p = 1:P
    
    % sample points
    n = [0:N_p(p), -N_p(p):-1].'; 
    
    % sampling matrix
    S_p = exp(1j*pi*n*k/K);
    S = blkdiag(S, S_p);
    
    % interpolation matrix
    Lambda_p = sinc(Delta_rel(:,p) - n.');
    Lambda = [Lambda, Lambda_p];
    
end

% total number of required samples
N = size(S,1);

% complexity of sampling (assuming use of iFFT)
C_samp = 8*P*K*log2(2*K);

% relative complexity
C_rel = (J*N + C_samp)/(2*J*PKmin1);

% approximation error
if return_eps_H
    H_approx = Lambda*S;
    eps_H = 10*log10(norm(H_approx - H, 'fro').^2/norm(H, 'fro').^2);    
end

end