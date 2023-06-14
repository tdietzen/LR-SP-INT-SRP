function [ Psi_STFT ] = calc_FD_GCC(y_STFT)
% [ Psi_STFT ] = calc_FD_GCC(y_STFT) computes the frequency domain GCC from
% STFT frames.
%
% IN:
% y_STFT        STFT tensor - freqbins x frames x microphone
%
% OUT:
% Psi_STFT      FD GCC tensor - freqbins x frames x microphone pairs

% get dimensions
[Kplus1,L,M] = size(y_STFT);
K = Kplus1-1;
P = M*(M-1)/2;

% init
Psi_STFT = zeros(K-1,L,P);

% compute FD GCC
for k = 2:K
    
    for l = 1:L
        
        y = squeeze(y_STFT(k,l,:));

        % FD GCC
        psi=zeros(P,1); 
        p = 0;
        for mprime = 1:M
            for m = mprime+1:M
                p = p+1;
                psi(p) = y(m)*y(mprime)';
            end
        end
        
       % PHAT weighting
        psi = psi./(abs(psi)+ 1e-9);
        Psi_STFT(k-1,l,:) = psi;
        
    end
end

end