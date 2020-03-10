%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2019 Thomas Dietzen
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt).
%
% If you find it useful, please cite:
%
% [1] T. Dietzen, S. Doclo, M. Moonen, and T. van Waterschoot, “Integrated
% sidelobe cancellation and linear prediction Kalman filter for joint 
% multimicrophone dereverberation, interfering speech cancellation, and 
% noise reduction,” IEEE/ACM Trans. Audio, Speech, Lang. Process., vol. 28,
% pp. 740 – 754, Jan. 2020.
% [2] T. Dietzen, S. Doclo, M. Moonen, and T. van Waterschoot, “Square 
% root-based multi-source early PSD estimation and recursive RETF update
% in reverberant environments by means of the orthogonal Procrustes
% problem,” IEEE/ACM Trans. Audio, Speech, Lang. Process., vol. 28,
% pp. 755 – 769, Jan. 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example of the ISCLP Kalman filter and its low-complexity variants as
% described in [1,2]. The code in contained main.m requests an SNR vaule by
% the user, loads microphone signals, estimates early PSDs and updates
% RETFs using [3] (see also https://github.com/tdietzen/SQRT-PSD-RETF),
% runs the ISCLP Kalman filter, and plays back the enhanced signal.


%% PREAMBLE

clear;
cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));
set(0,'DefaultFigureWindowStyle','docked');


%% CONFIGURATION

%%% ACOUSTIC SETTINGS

% speed of sound
c = 340;
% microphone positions
micPos = [...
    0, 0;...
    0.08, 0;...
    0.16, 0;...
    0.24, 0;...
    0.32, 0;...
    ];
% number of microphones
M = size(micPos,1);
% source angles
sourceAng = [0 60];
% SNR
SNR = input('which SNR/dB? ');
% microphone signal components
x1_TD           = audioread(['.' filesep 'audio' filesep 'x1.wav']);
s1_TD           = audioread(['.' filesep 'audio' filesep 's1.wav']);
x2_TD           = audioread(['.' filesep 'audio' filesep 'x2.wav']);
v_TD_SNR_0dB    = audioread(['.' filesep 'audio' filesep 'v_SNR_0dB.wav']);
v_TD_SNR_scaled = db2mag(-SNR)*v_TD_SNR_0dB;
y_TD  = x1_TD + x2_TD + v_TD_SNR_scaled;

%%% ALGORITHMIC SETTINGS

%%% STFT
% sample rate
fs = 16000;
% STFT parameters
N_STFT = 512;
R_STFT = N_STFT/2;
win = sqrt(hann(N_STFT,'periodic'));
N_STFT_half = floor(N_STFT/2)+1;
% frequency vector
f = linspace(0,fs/2,N_STFT_half);

%%% ISCLP Kalman filter [1,2]
% prediction length
L = 6;
% forgetting factor alpha
alpha_ISCLP_KF = 1-db2pow(-25);
A = sqrt(alpha_ISCLP_KF);
% LP filter variance
psi_wLP = db2pow(-4);
% SC filter variance
psi_wSC = db2pow(linspace(0,-15,257));
% build Psi_w_tilde_init and Psi_w_delta
psi_LP_init = kron((psi_wLP.^(1:L-1)).', ones(M,1));
Psi_w_tilde_init = cell(N_STFT_half,1);
Psi_w_delta = cell(N_STFT_half,1);
for k = 1:N_STFT_half
    Psi_w_tilde_init{k} = diag([psi_wSC(k)*ones(M-1,1); psi_LP_init]);
    Psi_w_delta{k} = (1-alpha_ISCLP_KF)*Psi_w_tilde_init{k};
end
% gain decay limitation
beta_SCLP_KF = db2mag(-2);

%%% PSD estimation and RETF update [3]
% forgetting factor zeta
zeta  = tau2forget(2*M*R_STFT/fs, R_STFT, fs);
% laplace coefficients of speech per frequency bin
tmp      = load('lap_div.mat');
lap_div  = tmp.lap_div;
% RETF updating threshold
xi_thresh = db2pow(-2);
% penelty factor alpha
tmp      = load('alpha_sqrtMP.mat');
alpha_SQRT_PSD_RETF = tmp.alpha_sqrtMP;
% initial RETFs
H_init_FT = doa2steervec(micPos, sourceAng, N_STFT_half, fs, c);
% diffuse coherence matrix
Gamma = calc_diffcoherence(micPos,N_STFT,fs,c,1e-3);


%%% FIGURE SETTINGS

% spectogram figure settings
xTickProp = [0, R_STFT/fs, fs/R_STFT];
yTickProp = [0, fs/(2000*R_STFT), R_STFT/2];
cRange    = [-55 5];


%% STFT PROCESSING

% transform
y_STFT  = calc_STFT(y_TD,  fs, win, N_STFT, R_STFT, 'onesided');
s1_STFT = calc_STFT(s1_TD, fs, win, N_STFT, R_STFT, 'onesided');

% plot
figure('Name',['microphone signal, SNR = ' num2str(SNR) 'dB']);
plotSpec(y_STFT(:,:,1),  'mag', xTickProp, yTickProp, cRange, 0); title(['y, SNR = ' num2str(SNR) 'dB']); ylabel('f/kHz');
figure('Name','target signal');
plotSpec(s1_STFT(:,:,1), 'mag', xTickProp, yTickProp, cRange, 0); title('s1'); ylabel('f/kHz');
drawnow;


%% EARLY PSD ESTIMATION, RETF UPDATE [3]

fprintf(' * estimate early PSDs, update RETFs [3]...\n');

% correlation matrix of microphone signal
Psi_y_STFT = estim_corrmat(y_STFT, zeta);
% compute GEVD
[P_STFT, lambda_STFT] = desmooth_GEVD(Psi_y_STFT, Gamma,...
    'lambdaMin', 0,...
    'forgetPSD', zeta);

[phi_s_hat,...
    phi_xl_hat,...
    H_hat_prior_STFT,...
    H_hat_post_STFT,...
    H_update_pattern]...
    = estim_PSD_RETF(P_STFT, lambda_STFT, Gamma, H_init_FT,...
    'method', 'square-root MP',...
    'itmax', 20,...
    'alpha', alpha_SQRT_PSD_RETF,...
    'beta', 20*lap_div.^2,...
    'xiThresh', xi_thresh...
    );

% plot
figure('Name',['target PSD estimate, SNR = ' num2str(SNR) 'dB']);
plotSpec(phi_s_hat(:,:,1), 'pow', xTickProp, yTickProp, cRange, 0); title(['phi s1 hat, SNR = ' num2str(SNR) 'dB']); ylabel('f/kHz');
drawnow;


%% ISCLP KALMAN FILTER [1,2]

% number of frames
numFrames = size(y_STFT,2);

% init outputs
q_STFT                   = cell(1,4);
e_prio_STFT              = cell(1,4);
e_post_STFT              = cell(1,4);
e_post_smooth_STFT       = cell(1,4);
q_STFT(1,:)              = {zeros(N_STFT_half, numFrames)};
e_prio_STFT(1,:)         = {zeros(N_STFT_half, numFrames)};
e_post_STFT(1,:)         = {zeros(N_STFT_half, numFrames)};
e_post_smooth_STFT(1,:)  = {zeros(N_STFT_half, numFrames)};

for variantIdx = 1:4
    
    switch variantIdx
        case 1
            complexity = 'L2M2';
        case 2
            complexity = 'LM2';
        case 3
            complexity = 'L2M';
        case 4
            complexity = 'LM';
    end
    fprintf([' * run O(' complexity ')-cost ISCLP Kalman filer [1,2]...\n']);

    for k = 2:N_STFT_half
        
        % reorganize data
        y_stack      = shiftdim(squeeze(y_STFT(k,:,:)),1);                % from (1 x numFrames x M)       to (M x numFrames)
        h_stack      = shiftdim(squeeze(H_hat_post_STFT(k,:,:,1)),1);     % from (1 x numFrames x M x 1)   to (M x numFrames)
        psi_sT_stack = shiftdim(squeeze(phi_s_hat(k,:,1)),1);             % from (1 x numFrames x 1)       to (1 x numFrames)
        
        % run ISCLP Kalman filter and spectral post processor
        [ q_stack, ...
            e_prio_stack,...
            e_post_stack,...
            e_post_smooth_stack ] = ...
            ISCLP(...
            A, Psi_w_delta{k},...
            y_stack,...
            psi_sT_stack,...
            h_stack,...
            Psi_w_tilde_init{k},...
            beta_SCLP_KF,...
            complexity);
        
        % save output
        q_STFT{variantIdx}(k,:)               = q_stack;
        e_prio_STFT{variantIdx}(k,:)          = e_prio_stack;
        e_post_STFT{variantIdx}(k,:)          = e_post_stack;
        e_post_smooth_STFT{variantIdx}(k,:)   = e_post_smooth_stack;
        
    end

    figure('Name',['enhanced signal, O(' complexity ')-cost, SNR = ' num2str(SNR) 'dB']);
    plotSpec(e_post_smooth_STFT{variantIdx}, 'mag', xTickProp, yTickProp, cRange, 0); title(['e post, O(' complexity ')-cost, SNR = ' num2str(SNR) 'dB']); ylabel('f/kHz');
    drawnow;
end

%% ISTFT PROCESSING

e_post_smooth_TD = cell(1,4);
for variantIdx = 1:4
    e_post_smooth_TD{variantIdx} = calc_ISTFT(e_post_smooth_STFT{variantIdx}, win, N_STFT, R_STFT, 'onesided');
end

%% WRITE AUDIO AND PLAY BACK

audiowrite(['.' filesep 'audio' filesep 'v_SNR_' num2str(SNR) 'dB.wav'], v_TD_SNR_scaled, fs);
audiowrite(['.' filesep 'audio' filesep 'y_SNR_' num2str(SNR) 'dB.wav'], y_TD, fs);
audiowrite(['.' filesep 'audio' filesep 'e_post_smooth_O(L2M2)_SNR_' num2str(SNR) 'dB.wav'],e_post_smooth_TD{1},fs);
audiowrite(['.' filesep 'audio' filesep 'e_post_smooth_O(LM2)_SNR_' num2str(SNR) 'dB.wav'],e_post_smooth_TD{2},fs);
audiowrite(['.' filesep 'audio' filesep 'e_post_smooth_O(L2M)_SNR_' num2str(SNR) 'dB.wav'],e_post_smooth_TD{3},fs);
audiowrite(['.' filesep 'audio' filesep 'e_post_smooth_O(LM)_SNR_' num2str(SNR) 'dB.wav'],e_post_smooth_TD{4},fs);

fprintf(' * play microphone signal...\n'); 
sound(y_TD(:,1), fs); pause(length(y_TD(:,1))/fs);
fprintf(' * play enhanced signal, O(LM)-cost...\n'); 
sound(e_post_smooth_TD{4}(:,1), fs);

fprintf('\nDONE.\n'); 
