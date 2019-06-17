function [phi_s_hat, eps_phi_s_rel] = solve_convMP_simple(H_hat, Psi_xe, phiMin, phi_s)
% function [phi_s_hat, eps_phi_s_rel] = solve_convMP_simple(H_hat, Psi_xe, phiMin, phi_s)
% solves simple conventional MP.
%
% IN:
% H_hat         RETFs - channels x sources
% Psi_xe        correlation matrix - channels x channels
% phiMin        PSD threshold
% phi_s         ground truth PSD estimates (optional) - sources x 1
%
% OUT:
% phi_s_hat     PSD estimate - sources x 1
% eps_phi_s_rel estimation error

N = size(H_hat, 2);

% cost function terms
A = abs(H_hat'*H_hat).^2;
b = -real(diag(H_hat'*Psi_xe*H_hat));

R = 1e-8*(trace(A)/N)*eye(N);
phi_s_hat  = max(-(A + R)\b, phiMin);

if nargin < 4
    eps_phi_s_rel = [];
else
    % error
    e_phi_s = sqrt(phi_s) - sqrt(phi_s_hat);
    eps_phi_s_rel = e_phi_s'*e_phi_s/sum(phi_s);
end

end

