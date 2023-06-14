%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 Thomas Dietzen
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt).
%
% If you use the code, please cite:
%
% @article{dietzen2023low,
%   author = {Dietzen, T. and De Sena, E. and van Waterschoot, T.},
%   title = {Low-Complexity Steered Response Power Mapping based on
%   Low-Rank and Sparse Interpolation},
%   journal = {tba},
%   volume = {tba},
%   number = {tba},
%   month = {tba},
%   year = {2023},
%   pages = {tba}
% }
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example of the SLRI SRP and SSPI SRP algorithm as described in [1]. The
% code contained in main.m computes the SRP matrix and its approximations
% for a NF of FF scenario, generates microphone signals, computes SRP maps
% and error metrics.
%
% Within this repository, the geometry of the microphone array, the search
% grid and the source locations are pre-defined. The diffuse noise signals
% are generated based on the microphone array geometry. If you would like
% to use different geometries, please provide corresponding RIRs and noise.
%
% In order to run this code, Matlab requires several GB of RAM (around 6GB
% in the NF scenario and 2GB in the FF scenario).


clear;

cd(fileparts(which(mfilename)));
addpath(genpath(pwd));



%% SET CONFIGURATION PARAMETERS

% acoustic parameters
scenario = 'NF';    % scenario (possible options are NF or FF)
c = 340;            % speed of sound
SNR = 0;            % SNR in dB
RT = 200;           % reverb time in ms (possible options are 0, 200, 400, 800)

% algorithmic parameters
C_rel_target = 10.^(-3.5:.2:-0.3);      % targeted relative complexities, cf. (61), corresponds to implemented C_rel up to round-off errors
N_aux = 2;                              % auxilary samples, cf. (31)

% algorithms to evaluate
alg_type = {'conv', 'LR', 'SI', 'SLRI', 'SSPI'};
alg_linestyles = {'--k', 'r', '<b', '-sb', '-ob'};  % for plotting

% evaluation parameters
numSrc = 32;                        % number of sources (must be <= 512)
numFrames = 16;                     % number of frames
numComp = length(C_rel_target);     % number of complexities
numAlg = length(alg_type);          % number of algorithms to evaluate
sigma = 6;                          % sigma in localization accuracy, cf. (65)
switch scenario
    case 'NF'
        eps_th = .2;                % eps_th in localization accuracy, cf. (65), in meter
    case 'FF'
        eps_th = 2.5;               % eps_th in localization accuracy, cf. (65), in degree    
end



%% LOAD PRE-DEFINED CONFIGURATION PARAMETERS

% geometric parameters
geomParam = load(['./data/param/geomParam_' scenario '.mat']);
 
roomDim = geomParam.roomDim;    % room dimensions
micPos = geomParam.micPos;      % micorphone locations

M = size(micPos,1);             % number of mics
P = M*(M-1)/2;                  % number of mic pairs, cf. (6)

switch scenario
    
    case 'NF'
        coord_x = geomParam.coord_x;    % x-coordinates of candidate locations
        coord_y = geomParam.coord_y;    % y-coordinates of candidate locations
        coord_z = geomParam.coord_z;    % z-coordinates of candidate locations
        
    case 'FF'
        ang_az = geomParam.ang_az;  	% azimuth angles of candidate locations
        ang_pol = geomParam.ang_pol;    % polar angles of candidate locations
end

% freuqncy parameters
freqParam = load(['./data/param/freqParam_' scenario '.mat']);

fs = freqParam.fs; % sampling frequency
K = freqParam.K;   % frequency bins (one-sided)

win = sqrt(hann(2*K,'periodic'));   % window for STFT processing
frameOffset = ceil((.8*fs)/K);      % frame offset (number of frames corresponging to the length of the RIR, first and last frameOffset frames will not be used for evaluation)

srcSig = audioread(['./data/sig/srcSig_' scenario '.wav']);       % source signal
noiseSig = audioread(['./data/sig/noiseSig_' scenario '.wav']);   % noise signal

RIR_str = ['./data/RIR/RIR_' scenario '_'];                       % sting for RIRs



%% INIT PERFORMANCE MEASURES

% init performance measures
C_rel = cell(numAlg,1); % relative complexity, corresponds to C_rel_target up to round-off errors
eps_H = cell(numAlg,1); % matrix error
eps_z = cell(numAlg,1); % map error
eps_s = cell(numAlg,1); % localization error
rho_s = cell(numAlg,1); % localization accuracy

for idx = 1:numAlg
    switch alg_type{idx}
        case {'conv', 'SI'}
            C_rel{idx} = nan;
            eps_H{idx} = nan;
            eps_z{idx} = nan(numSrc, numFrames);
            eps_s{idx} = nan(numSrc, numFrames);
            rho_s{idx} = nan(numSrc, numFrames);
        case {'LR', 'SLRI', 'SSPI'}
            C_rel{idx} = nan(numComp, 1);
            eps_H{idx} = nan(numComp, 1);
            eps_z{idx} = nan(numSrc, numFrames, numComp);
            eps_s{idx} = nan(numSrc, numFrames, numComp);
            rho_s{idx} = nan(numSrc, numFrames, numComp);
    end
end



%% COMPUTE CANDIDATE LOCATIONS AND TDOAs


% get candidate locations q_i and TDOAs Delta_t_p
switch scenario
    case 'NF'
        [q_i, Delta_t_p] = gen_candloc(micPos, c, coord_x, coord_y, coord_z);
    case 'FF'
        [q_i, Delta_t_p] = gen_candloc(micPos, c, ang_pol, ang_az);
end


%%% plot geometry %%%

figure;
hold on;

% plot microphone locations
plot3(micPos(:,1), micPos(:,2), micPos(:,3), 'or');

% plot candidate locations
arrayCenter = mean(micPos,1);
switch scenario
    case 'NF'
        plot3(q_i(:,1), q_i(:,2), q_i(:,3), '.b');
    case 'FF'
        plot3(-2*q_i(:,1) + arrayCenter(1), -2*q_i(:,2) + arrayCenter(2), -2*q_i(:,3) + arrayCenter(3), '.');
        % since cand_loc corresponds to an incident direction vector
        % pointing towards the array, we apply a minus in order to obtain a
        % vector pointing out of the array.
end

% axis handling
axis equal
xlim([0 roomDim(1)]);
ylim([0 roomDim(2)]);
zlim([0 roomDim(3)]);
grid;

% label handling
xlabel('depth (m)');
ylabel('width (m)');
zlabel('height (m)');
legend('microphones', 'candidate locations');
title(['microphone and candidate locations in the ' scenario ' scenario']);

% adjust view angle
view([-40 30]);

drawnow;



%% COMPUTE MATRICES

w_b = pi*fs;         % bandlimit
w = w_b*(1:K-1).'/K; % omega


fprintf('\n *** BUILD SRP MATRICES *** \n\n');

%%% conventional SRP %%%

fprintf('Conventional SRP matrix H... '); tic;

matr_data{1}.H = exp(1j*kron(Delta_t_p, w.')); % cf. (19), (20)

fprintf(['took ' num2str(round(toc, 1)) 's.\n']);


%%% SRP approximations %%%

[ T, N_p ] = sample_param(micPos, w_b, c);   % sampling period and required samples required in SI, SLRI, and SSPI approximation
N_p = N_p + N_aux;                           % add auxiliary samples 

for idx = 2:numAlg
    
    fprintf([alg_type{idx} ' approximation of H... ']); tic;

    switch alg_type{idx}
        
        
        case 'LR'
            [...
                matr_data{idx}.H_tall,...
                matr_data{idx}.H_fat,...
                C_rel{idx},...
                eps_H{idx}] = lr_approx(matr_data{1}.H, C_rel_target); % cf. (25), (26)
            
        case 'SI'
            [...
                matr_data{idx}.Lambda,...
                matr_data{idx}.S,...
                C_rel{idx},...
                eps_H{idx}] = si_approx(matr_data{1}.H, Delta_t_p/T, N_p); % cf. (35), (37), (40), (41)
            
        case 'SLRI'
            [...
                matr_data{idx}.Lambda_tall,...
                matr_data{idx}.Lambda_fat,...
                ~,...
                C_rel{idx},...
                eps_H{idx}] = slci_approx(matr_data{1}.H, Delta_t_p/T, N_p, C_rel_target, 'LR'); % cf. (52)
            
        case 'SSPI'
            [...
                matr_data{idx}.Lambda_sp,...
                ~,...
                C_rel{idx},...
                eps_H{idx}] = slci_approx(matr_data{1}.H, Delta_t_p/T, N_p, C_rel_target, 'SP'); % cf. (53)
            
    end
    fprintf(['took ' num2str(round(toc, 1)) 's.\n']);

end

%%% plot eps_H %%%

figure;
hold on;

% plot errors
for idx = 2:numAlg
    plot(C_rel{idx}, eps_H{idx}, alg_linestyles{idx});
end

% axis handling
set(gca, 'Xscale', 'log');
xlim([min(C_rel_target) max(C_rel_target)]);
ylim([-35 0]);
grid;

% label handling
xlabel('relative complexity');
ylabel('matrix error (dB)');
legend(alg_type(2:end));
title(['matrix error in the ' scenario ' scenario']);

drawnow;



%% COMPUTE MAPS

fprintf(['\n *** COMPUTE SRP MAPS *** \n\n']);

fprintf(['SNR = ' num2str(SNR) 'dB \n']);
fprintf(['RT = ' num2str(RT) 'ms \n\n']);

for srcIdx = 1:numSrc
    
    
    %%% generate data %%%
    
    % load RIR data
    rir_data = load([RIR_str 'src' num2str(srcIdx) '_rt' num2str(RT) '.mat']);
    rir_dry_data = load([RIR_str 'src' num2str(srcIdx) '_rt0.mat']); % dry signals needed as a reference for SNR
    
    % define true source location
    switch scenario
        case 'FF'
            q_s = mean(micPos,1) - rir_data.q_s; % incident vector pointing towards array
            q_s = q_s/norm(q_s);
        case 'NF'
            q_s = rir_data.q_s;
    end
    
    % convolve source signal with RIRs
    x = zeros(size(rir_data.RIR, 1) + length(srcSig) - 1, M);           % reverberant signal
    x_dry = zeros(size(rir_dry_data.RIR, 1) + length(srcSig) - 1, M);   % dry signals needed as a reference for SNR
    for m = 1:M
        x(:,m) = fast_conv(rir_data.RIR(:,m), srcSig);
        x_dry(:,m) = fast_conv(rir_dry_data.RIR(:,m), srcSig);
    end
    
    % scale noise signal
    v = noiseSig(1:size(x,1),:);
    scaling = scaling_SNR(x_dry, v, SNR);
    v = scaling*v;
    
    % transform to STFT domain
    x_STFT = calc_STFT(x, fs, win, 2*K, K, 'onesided');
    v_STFT = calc_STFT(v, fs, win, 2*K, K, 'onesided');
    
    % select frames with sufficient signal level
    rng(srcIdx);
    select_thresh = SNR - 3;
    y_STFT = select_frames(x_STFT, v_STFT, numFrames, frameOffset, select_thresh);
    
    % compute FD GCC
    psi_STFT = calc_FD_GCC(y_STFT);   % cf. (2)
    
    
    fprintf(['source ' num2str(srcIdx) '... ']);
    tic;
    
    for frameIdx = 1:numFrames
        
        % get FD GCC vector psi and TD GCC sample vector xi
        psi = [];
        xi = [];
        for p = 1:P
            
            % load FD GCC
            psi_p = psi_STFT(:, frameIdx, p);
            
            % compute TD GCC samples using iFFT
            psi_p_bar = [0; psi_p; 0; conj(psi_p(end:-1:1))];   % cf. (67)
            n = [0:N_p(p), -N_p(p):-1].';
            n_prime_mod =(mod(n,2*K)+1).';                      % cf. (46)
            xi_p_bar = ifft(psi_p_bar, 'symmetric');            % cf. (44)
            xi_p = 2*K*xi_p_bar(n_prime_mod);                   % cf. (45)
            
            % stack FD GCCs
            psi = [psi; psi_p];    % cf. (20)
            
            % stack TD GCC samples
            xi = [xi; xi_p];       % cf. (36)
            
        end
        
        %%% compute maps %%%
        
        for idx = 1:numAlg
            
            switch alg_type{idx}
                
                case 'conv'
                    
                    % load matrix
                    H = matr_data{idx}.H;
                    % compute map
                    z = 2*real(H*psi);                                                                  % cf. (22)
                    % compute localization error
                    [eps_s{idx}(srcIdx, frameIdx), src_loc_hat] = loc_err(z, q_i, q_s, scenario);   % cf. (64)
                
                case 'LR'
                    
                    % load matrices
                    H_tall = matr_data{idx}.H_tall;
                    H_fat = matr_data{idx}.H_fat;
                    
                    for compIdx = 1:numComp
                        
                        % compute map
                        z_type = 2*real(H_tall{compIdx}*(H_fat{compIdx}*psi));                              % cf. (27)
                        % compute map error
                        eps_z{idx}(srcIdx, frameIdx, compIdx) = 10*log10(norm(z_type - z).^2/norm(z).^2);   % cf. (63)
                        % compute localization error
                        eps_s{idx}(srcIdx, frameIdx, compIdx) = loc_err(z_type, q_i, q_s, scenario);    % cf. (64)
                        
                    end
                
                case 'SI'
                    
                    % load matrix
                    Lambda = matr_data{idx}.Lambda;
                    % compute map
                    z_type = Lambda*xi;                                                         % cf. (42)
                    % compute map error
                    eps_z{idx}(srcIdx, frameIdx) = 10*log10(norm(z_type - z).^2/norm(z).^2);    % cf. (63)
                    % compute localization error
                    eps_s{idx}(srcIdx, frameIdx) = loc_err(z_type, q_i, q_s, scenario);     % cf. (64)
                
                case 'SLRI'
                    
                    % load matrices
                    Lambda_tall = matr_data{idx}.Lambda_tall;
                    Lambda_fat = matr_data{idx}.Lambda_fat;
                    
                    for compIdx = 1:numComp
                        if ~isempty(Lambda_tall{compIdx})
                            % compute map
                            z_type = Lambda_tall{compIdx}*(Lambda_fat{compIdx}*xi);                             % cf. (55)
                            % compute map error
                            eps_z{idx}(srcIdx, frameIdx, compIdx) = 10*log10(norm(z_type - z).^2/norm(z).^2);   % cf. (63)
                            % compute localization error
                            eps_s{idx}(srcIdx, frameIdx, compIdx) = loc_err(z_type, q_i, q_s, scenario);    % cf. (64)
                        end
                    end
                
                case 'SSPI'
                    
                    % load matrix
                    Lambda_sp = matr_data{idx}.Lambda_sp;
                    
                    for compIdx = 1:numComp
                        if ~isempty(Lambda_sp{compIdx})
                            % compute map
                            z_type = Lambda_sp{compIdx}*xi;                                                     % cf. (56)
                            % compute map error
                            eps_z{idx}(srcIdx, frameIdx, compIdx) = 10*log10(norm(z_type - z).^2/norm(z).^2);   % cf. (63)
                            % compute localization error
                            eps_s{idx}(srcIdx, frameIdx, compIdx) = loc_err(z_type, q_i, q_s, scenario);    % cf. (64)
                        end
                    end
                    
            end
        end
    end
    fprintf(['took ' num2str(round(toc, 1)) 's.\n']);
    
end


%%% plot eps_z %%%

figure;
hold on;

% plot map errors
for idx = 2:numAlg
    plot(C_rel{idx}, median(reshape(eps_z{idx}, [], length(C_rel{idx}))), alg_linestyles{idx});
end


% axis handling
set(gca, 'Xscale', 'log');
xlim([min(C_rel_target) max(C_rel_target)]);
ylim([-35 0]);
grid;

% label handling
xlabel('relative complexity');
ylabel('median map error (dB)');
legend(alg_type(2:end));
title(['map error in the ' scenario ' scenario']);

drawnow;


%%% plot rho_s %%%

figure;
hold on;

% plot accuracy
for idx = 1:numAlg
    
    % convert to localization accuraracy
    rho_s{idx} = soft_detect(eps_s{idx}, sigma, eps_th);
    
    if idx == 1
        plot([min(C_rel_target) max(C_rel_target)], [1 1]*mean(rho_s{idx}(:)), alg_linestyles{idx});
    else
        plot(C_rel{idx}, mean(reshape(rho_s{idx}, [], length(C_rel{idx}))), alg_linestyles{idx});
    end
end

% axis handling
set(gca, 'Xscale', 'log');
xlim([min(C_rel_target) max(C_rel_target)]);
ylim([0 1]);
grid;

% label handling
xlabel('relative complexity');
ylabel('localization accuracy');
legend(alg_type(1:end));
title(['localization accuracy in the' scenario ' scenario']);

drawnow;

fprintf(['\n DONE. \n\n']);





%% HELPER FUNCTIONS

function [y_STFT] = select_frames(x_STFT, v_STFT, numFrames, framesOffset, select_thresh)

% select frames with sufficient signal level
L = size(x_STFT,2);
useframe_idx = false(L,1);
for l = framesOffset+1:size(x_STFT,2)-framesOffset
    SNR_local = mag2db(norm(squeeze(x_STFT(:,l,:)), 'fro')/(norm(squeeze(v_STFT(:,l,:)), 'fro') + eps));
    % select frames only with sufficient signal level
    if SNR_local > select_thresh
        useframe_idx(l) = true;
    end
end
% convert logical index to integer index
useframe_idx = find(useframe_idx);
usableFrames = length(useframe_idx);
% check if we have enough frames
if usableFrames < numFrames
    error('too few frames selected.');
end
% randomly select
useframe_idx = useframe_idx(randperm(usableFrames,numFrames));

y_STFT = x_STFT(:,useframe_idx,:) + v_STFT(:,useframe_idx,:);

end