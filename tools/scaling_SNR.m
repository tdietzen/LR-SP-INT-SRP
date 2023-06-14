function [scaling] = scaling_SNR(x, v, SNR)
% [scaling] = scaling_SNR(x, v, SNR)
% get scaling for noise v relative to signal x to obtain specified SNR.
%
% IN:
% x           speech signal
% v           noise signal
% SNR         desired SNR
%
% OUT:
% scaling     scaling factor


% get powers
power_x = norm(x, 'fro');
power_v  = norm(v, 'fro');

% get scaling
if isinf(SNR)
    scaling = 0;
else
    scaling = ((power_x/power_v)/db2mag(SNR));
end

end