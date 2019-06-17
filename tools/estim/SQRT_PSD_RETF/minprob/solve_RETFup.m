function [H_hat_post, up, rho_H_prior_rel, rho_H_post_rel] = solve_RETFup( sqrtPsi_xe, sqrtphi_s_hat, Omega_hat, xiThresh, reg, beta, H_hat_prior, H)
% [H_hat_post, up, rho_H_prior_rel, rho_H_post_rel] = solve_RETFup( sqrtPsi_xe, sqrtphi_s_hat, Omega_hat, xiThresh, reg, beta, H_hat_prior, H)
% solves RETF update MP.
%
% IN:
% sqrtPsi_xe        square root of correlation matrix - channels x sources
% sqrtphi_s_hat     square root PSD estimate - sources x 1
% Omega_hat         estimate of unitary matrix - sources x sources
% xiThresh          update threshold
% reg               threshold regularization
% beta              penalty factor 
% H_hat_prior       prior RETFs - channels x sources
% H                 ground truth RETF (optional)
%
% OUT:
% H_hat_post        posterior RETFs - channels x sources
% up                RETF update pattern - channels x sources
% rho_H_prior_rel   prior estimation error
% rho_H_post_rel    posterior estimation error

[M, N] = size(H_hat_prior);

% init
H_hat_post = ones(M,N);

% update
phi_s_hat = real(diag(sqrtphi_s_hat')*sqrtphi_s_hat);
xi = phi_s_hat/(sum(phi_s_hat) + reg);

up = zeros(N,1);
for n = 1:N
    if xi(n) > xiThresh
        up(n) = 1;
        H_hat_post(2:M,n) = (sqrtPsi_xe(2:M,:)*Omega_hat(:,n)*conj(sqrtphi_s_hat(n)) + beta*H_hat_prior(2:M,n))/(phi_s_hat(n) + beta);
    else
        up(n) = 0;
        H_hat_post(2:M,n) = H_hat_prior(2:M,n);
    end
end





%if nargin < 6
if nargin < 8
    rho_H_prior_rel = [];
    rho_H_post_rel  = [];
else 
    % relative error
    H_e_prior   = H - H_hat_prior;
    H_e_post    = H - H_hat_post;
    
    rho_H_prior_rel = trace(H_e_prior'*H_e_prior)/trace(H(2:M,:)'*H(2:M,:));
    rho_H_post_rel  = trace(H_e_post'*H_e_post)/trace(H(2:M,:)'*H(2:M,:));
end

end          
