function [y] = soft_detect(x, sigma, eps_th)
% [y] = soft_detect(x, sigma, eps_th) maps x between 0 and 1.
%
% IN:
% x       data to be mapped
% sigma   hardness of the threshold
% eps_th  threshold
%
% OUT:
% y       mapped data

y = 1./(1+exp((sigma/eps_th)*(x-eps_th)));