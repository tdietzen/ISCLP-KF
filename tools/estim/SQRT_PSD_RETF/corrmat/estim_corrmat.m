function [Psi_x_smth, Psi_x_mean] = estim_corrmat(X, alpha)
% [Psi_x_smth, Psi_x_mean] = estim_corrmat(X, alpha) 
% estimates correlation matrix.
%
% IN:
% X             STFT data - freqbins x frames x channels
% alpha         forgetting factor
% 
% OUT
% Psi_x_smth    smooth correlation matrix estimate - freqbins x frames x channels x channels
% Psi_x_mean    mean correlation matrix estimate - freqbins x 1 x channels x channels

% dimensions
[N_half, L, M, ~] = size(X);          

%%% compute average PSD matrix
if alpha ~= 1
    R_inst = zeros(N_half, L, M, M);
    for l = 1:L
        for k = 1:N_half
            x = squeeze(X(k,l,:));
            R_inst(k,l,:,:) = shiftdim(x*x', -2);
        end
    end
    Psi_x_mean = sum(R_inst, 2)/L;
else
    Psi_x_mean = zeros(N_half, 1, M, M);
    for l = 1:L
        for k = 1:N_half
            x = squeeze(X(k,l,:));
            Psi_x_mean(k,1,:,:) = Psi_x_mean(k,1,:,:) + shiftdim(x*x', -2);
        end
    end
    Psi_x_mean = Psi_x_mean/L;
end

%%% compute smooth PSD matrix

R_smth_tmp = eps*ones(N_half, 1, M, M);
if alpha ~= 1
    Psi_x_smth = zeros(N_half, L, M, M); 
    for l = 1:L
        R_smth_tmp = alpha*R_smth_tmp + (1-alpha)*R_inst(:,l,:,:);
        Psi_x_smth(:,l,:,:) = R_smth_tmp;
    end
else
    Psi_x_smth = Psi_x_mean;
end
