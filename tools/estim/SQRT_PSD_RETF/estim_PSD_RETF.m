function [phi_s_hat_STFT, phi_d_STFT, H_hat_prior_STFT, H_hat_post_STFT, update_RETF] = estim_PSD_RETF(P_STFT, lambda_STFT, Gamma_FT, H_init_FT, varargin) 
% function [phi_s_hat_STFT, phi_d_STFT, H_hat_prior_STFT, H_hat_post_STFT, update_RETF] = estim_PSD_RETF(P_STFT, lambda_STFT, Gamma_FT, H_init_FT, varargin) 
% estimates PSDs and updates RETFs.
%
% IN:
% P_STFT                    eigenvectors - freqbins x frames x channels x channels
% lambda_STFT               eigenvalues - freqbins x frames x channels
% Gamma_FT                  diffuse coherence matrix - freqbins x 1 x frames x channels
% H_init_FT                 initial RETF estimate - freqbins x frames x channels x sources
% 'itmax', itmax            maximum number of iterations
% 'phiMin', phiMin          power threshold
% 'method', method          {'conventional MP', 'square-root MP'}, defines estimation method
% 'alpha', alpha            penalty factor conventional and square-root MP
% 'beta', beta              penalty factor in RETF update MP
% 'xiThresh', xi_thresh     threshold for RETF update
%
% OUT:
% phi_s_hat_STFT            early PSD estimates - freqbins x frames x sources
% phi_d_STFT                diffuse PSD estimate - freqbins x frames
% H_hat_prior_STFT          prior RETF estimate - freqbins x frames x channels x sources
% H_hat_post_STFT           posterior RETF estimate - freqbins x frames x channels x sources
% update_RETF               RETF update pattern - freqbins x frames x sources

% dimensions
[N_FT_half, L, M]      = size(lambda_STFT);             % number of frequency bins, frames, microphones
N                      = size(H_init_FT, 4);            % number of sources

% default values
method                    = 'square-root MP';
itmax                     = 20;                         % max iterations
phiMin                    = db2pow(-80);                % minimum power threshold
beta                      = 1e3*ones(N_FT_half,1);      % beta
alpha                     = 1e3*ones(N_FT_half,1);      % alpha

% parse name-value pairs
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        switch varargin{i} 
            case 'itmax'
                itmax       = varargin{i+1}; 
            case 'phiMin'
                phiMin      = varargin{i+1};          
            case 'method'    
                method      = varargin{i+1};   
            case 'alpha'
                alpha       = varargin{i+1};
            case 'beta'
                beta        = varargin{i+1};   
            case 'xiThresh'
                xi_thresh   = varargin{i+1};                   
        end
    end
end

% init
phi_s_hat_STFT        = zeros(N_FT_half,L,N);
phi_d_STFT             = zeros(N_FT_half,L);

if strcmp('method','square-root MP')
    H_hat_prior_STFT       = zeros(N_FT_half,L,M,N);
    H_hat_post_STFT        = zeros(N_FT_half,L,M,N);
    update_RETF            = zeros(N_FT_half,L,N);
else
    H_hat_prior_STFT = [];
    H_hat_post_STFT  = [];
    update_RETF      = [];
end


for k=2:N_FT_half
    
    % get H
    H_hat = squeeze(H_init_FT(k,1,:,:));
    H_hat_post = H_hat;
    
    % get Gamma
    Gamma = squeeze(Gamma_FT(k,1,:,:));
    
    for l = 1:L
               
        % get P and lambda
        P = squeeze(P_STFT(k,l,:,:));
        lambda = squeeze(lambda_STFT(k,l,:));
        
        %%% Decompose into early and late part %%%
        %
        % find index of maximum eigenvalues and principal eigenvectors
        [~, maxIdx] = sort(lambda, 'descend');
        Pmax = P(:,maxIdx(1:N)); 
        lambdaMax = lambda(maxIdx(1:N));
        lambdaMin = lambda(maxIdx(N+1:end));
        %
        % compute late reverberant PSD
        phi_d = mean(lambdaMin); % using min or mean doesn't make much of a difference
        %
        % save
        phi_d_STFT(k,l) = phi_d;
        %
        % compute lambda_xe
        lambda_xe = lambdaMax - phi_d;
        lambda_xe(lambda_xe < phiMin) = phiMin;

        
        %%% early correlation matrix %%%
        %
        sqrtlambda_xe = sqrt(lambda_xe);
        sqrtPsi_xe = Gamma*Pmax*diag(sqrtlambda_xe);
        Psi_xe = sqrtPsi_xe*sqrtPsi_xe';
        
        
        switch method

            case 'conventional MP'
                
                % solve
                phi_s_hat = solve_convMP(H_hat, Psi_xe, phiMin, alpha(k), itmax);

            case 'square-root MP'
                
                % propagate H_hat_post
                H_hat_prior = H_hat_post;
                
                % initial value
                sqrtphi_s_init = sqrt(solve_convMP_simple(H_hat_prior, Psi_xe, phiMin));
                
                % solve
                [sqrtphi_s_hat, Omega_hat] = solve_sqrtMP(sqrtPsi_xe, H_hat_prior, sqrtphi_s_init, alpha(k), itmax);
                phi_s_hat = sqrtphi_s_hat.*conj(sqrtphi_s_hat);
                
                %%% update RETF H %%%
                if l > 16
                    % update H
                    phi_reg = phi_d + 1e-3;
                    [H_hat_post, up] = solve_RETFup(sqrtPsi_xe, sqrtphi_s_hat, Omega_hat, xi_thresh, phi_reg, beta(k), H_hat_prior);

                else
                    H_hat_post = H_hat_prior;
                    up = zeros(N,1);
                end
                
                % save
                H_hat_prior_STFT(k,l,:,:) = shiftdim(H_hat_prior, -2);
                H_hat_post_STFT(k,l,:,:) = shiftdim(H_hat_post, -2);
                update_RETF(k,l,:) = up;
        end
        
        % save
        phi_s_hat_STFT(k,l,:) = phi_s_hat;

    end   
end


end