function [sqrtphi_s_hat, Omega_hat, eps_phi_s_rel, eps_phi_s_rel_it] = solve_sqrtMP(sqrtPsi_xe, H_hat, sqrtphi_s_init, alpha, itMax, phi_s)
% [sqrtphi_s_hat, Omega_hat, eps_phi_s_rel, eps_phi_s_rel_it] = solve_sqrtMP(sqrtPsi_xe, H_hat, sqrtphi_s_init, alpha, itMax, phi_s)
% solves square-root MP.
%
% IN:
% sqrtPsi_xe        square root of correlation matrix - channels x sources
% H_hat             RETFs - channels x sources
% sqrtphi_s_init    intital square root PSD estimate - sources x 1
% alpha             penalty factor
% itMax             max number of iterations
% phi_s             ground truth PSD estimates (optional) - sources x 1
%
% OUT:
% phi_s_hat         square root PSD estimate - sources x 1
% Omega_hat         estimate of unitary matrix - sources x sources
% eps_phi_s_rel     estimation error
% eps_phi_s_rel_it  estimation error per iteration

% initial value
sqrtphi_s_hat = sqrtphi_s_init;

% number of sources
[M, N] = size(H_hat);


if nargin < 6
    compute_eps = 0;
else
    compute_eps = 1;
end

if compute_eps
    eps_phi_s_rel_it = zeros(itMax, 1);
else
    eps_phi_s_rel_it = [];
    eps_phi_s_rel = [];
end

for i_it = 1:itMax
    
    % save convergence plot
    if compute_eps
        e_phi_s = sqrt(phi_s) - abs(sqrtphi_s_hat);
        eps_phi_s_rel_it(i_it) = e_phi_s'*e_phi_s/sum(phi_s);
    end
    
    % get Sigma
    [U, ~, V] = svd(sqrtPsi_xe'*H_hat*diag(sqrtphi_s_hat));

    Omega_hat = U*V';
    
    % get sqrtphi_xe_hat
    sqrtphi_s_hat_old = sqrtphi_s_hat;
    
    sqrtphi_s_hat = (diag(H_hat'*sqrtPsi_xe*Omega_hat) + alpha*Omega_hat'*sqrtPsi_xe(1,:)')./...
        (diag((H_hat'*H_hat)) + alpha*ones(N,1));
    
    % break condition
    delta_sqrtphi_s_hat = abs(sqrtphi_s_hat) - abs(sqrtphi_s_hat_old);
    if delta_sqrtphi_s_hat'*delta_sqrtphi_s_hat/(sqrtphi_s_hat_old'*sqrtphi_s_hat_old) < 1e-6
        % save convergence plot and break
        if compute_eps
            eps_phi_s_rel_it(i_it+1:end) = eps_phi_s_rel_it(i_it);
        end
        break;
    end
end

if compute_eps
    eps_phi_s_rel = eps_phi_s_rel_it(end);
end

end

