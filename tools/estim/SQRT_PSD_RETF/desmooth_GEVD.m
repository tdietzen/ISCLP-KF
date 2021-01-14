function [P_STFT, lambda_STFT, lambda_LP_STFT] = desmooth_GEVD(Psi_x_STFT, Gamma_FT, varargin) 
% [P_STFT, lambda_STFT, lambda_LP_STFT] = desmooth_GEVD(Psi_x_STFT, Gamma_FT, varargin) 
% performs GEVD and subspace-based desmoothing.
%
% IN:
% Psi_x_STFT                correlation matrix - freqbins x frames x channels x channels
% Gamma_FT                  diffuse coherence matrix - freqbins x 1 x channels x channels
% 'lambdaMin', lambdaMin    eigenvalue threshold
% 'forgetPSD', forgetPSD    forgetting factor for desmoothing
%
% OUT:
% P_STFT                    eingevectors - freqbins x frames x channels x channels
% lambda_STFT               desmoothed eigenvalues - freqbins x frames x channels
% lambda_LP_STFT            smooth eigenvalues - freqbins x frames x channels


% dimensions
[N_FT_half, L, M, ~]      = size(Psi_x_STFT);     % number of frequency bins, frames, microphones

% default options
forgetPSD  = 0;                          % PSD forgetting factor
lambdaMin  = 0;                          % minimum power

% read options from input
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}              
            case 'forgetPSD'
                forgetPSD     = varargin{i+1};              
            case 'lambdaMin'
                lambdaMin     = varargin{i+1};              
        end
    end
end

% init
P_STFT                 = zeros(N_FT_half,L,M,M);
lambda_STFT            = zeros(N_FT_half,L,M);
lambda_LP_STFT         = zeros(N_FT_half,L,M);

for k=2:N_FT_half
    R_v = 0;
    
    % apply regularization
    Gamma = squeeze(Gamma_FT(k,1,:,:));
    for l = 1:L
        
        %%% GEVD %%%
        %
        Psi_x = squeeze(Psi_x_STFT(k,l,:,:));
        % generalized eigenvalue decomposition
        [P, Lambda_LP] = eig(Psi_x, Gamma);
        % ignore residual imaginary component and set negative values to zero
        lambda_LP = real(diag(Lambda_LP));
        % rescaling W such that P'*Rv_k*P = I and P'*Ry_kp*P = Lambda
        P = P/diag(sqrt(real(diag(P'*Gamma*P))));
        
        
        %%% Sorting eigenvalues/eigenvectors %%%
        % 
        if l > 1
            % map current eigenvector order to previous one, start with principal eigenvector 
            [~, maxIdx] = sort(lambda_LP, 'descend');
            % if ordered correctly, Q approximates I
            Q = abs(P_old'*Gamma*P);
            % temporary variables
            P_tmp = P;
            lambda_LP_tmp = lambda_LP;
            % order such that Q approximates I
            for m = transpose(maxIdx)
                % find index
                [~, maxIdx] = max(Q(:,m));
                % set corresponding row to zero (guarantees that maxIdx is unique in each iteration) 
                Q(maxIdx,:) = 0;
                % sort P and Lambda
                P(:, maxIdx) = P_tmp(:,m);
                lambda_LP(maxIdx) = lambda_LP_tmp(m);
            end
        end
        P_old = P;    
        %
        % save 
        P_STFT(k,l,:,:)       = shiftdim(P, -2);
        lambda_LP_STFT(k,l,:) = shiftdim(lambda_LP, -2);
        
        
        %%% Filter eigenvalues to compensate for recursive avaeraging %%%
        %
        if l > 1
            % apply HP to Lambda
            lambda_LP_old = squeeze(lambda_LP_STFT(k,l-1,:));
            lambda = 1/(1-forgetPSD)*lambda_LP - forgetPSD/(1-forgetPSD)*lambda_LP_old;
            lambda(lambda < lambdaMin) = lambdaMin;
        else
            % use LP version
            lambda = lambda_LP;
            lambda(lambda < lambdaMin) = lambdaMin;
        end
        %
        % save
        lambda_STFT(k,l,:) = shiftdim(lambda, -2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % NOTE: instead of desmoothing eigenvalues as in the implementation
        % above, one might prefer to use instantaneous eigenvalues, see
        % also [a] and the desmooth_GEVD() function in [b].
        %
        % [a] T. Dietzen, M. Moonen, and T. van Waterschoot, 'Instantaneous
        % PSD estimation for speech enhancement based on generalized
        % principal components,' in Proc. 28th European Signal Process.
        % Conf. (EUSIPCO 2020), Amsterdam, Netherlands, Jan 2021, pp. 1-5.
        %
        % [b] https://github.com/tdietzen/INSTANT-PSD
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end               
end

end
