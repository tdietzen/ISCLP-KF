function [phi_s_hat, eps_phi_s_rel] = solve_convMP(H_hat, Psi_xe, phiMin, alpha, itMax, phi_s)
% function [phi_s_hat, eps_phi_s_rel] = solve_convMP(H_hat, Psi_xe, phiMin, alpha, itMax, phi_s)
% solves conventional MP.
%
% IN:
% H_hat         RETFs - channels x sources
% Psi_xe        correlation matrix - channels x channels
% phiMin        PSD threshold
% alpha         penalty factor
% itMax         max number of iterations
% phi_s         ground truth PSD estimates (optional) - sources x 1
%
% OUT:
% phi_s_hat     PSD estimate - sources x 1
% eps_phi_s_rel estimation error


N = size(H_hat,2);

% cost function terms
A1 = abs(H_hat'*H_hat).^2;
A2 = alpha*ones(N);
A = A1 + A2;
b = -real(diag(H_hat'*Psi_xe*H_hat) + alpha*Psi_xe(1,1)*ones(N,1));

% initial value
R = 1e-8*(trace(A1)/N)*eye(N);
x0  = -(A + R)\b;

% gradient
G   = @(x) (x'*A)'+b;

% step size
mu = 1/norm(A);

for it = 1:itMax
    
    xOld = x0;
    % compute gradient
    Gtemp = G(x0);
    % take gradient step
    xstep = x0-mu*Gtemp;
    % apply proximal opetor
    xProj = max(xstep,phiMin);
    % overwrite
    x0 = xProj;
    
    if norm(x0-xOld)/norm(xOld)<1e-6
        break;
    end
    
end

phi_s_hat = x0;


if nargin < 6
    eps_phi_s_rel = [];
else
    % error
    e_phi_s = sqrt(phi_s) - sqrt(phi_s_hat);
    eps_phi_s_rel = e_phi_s'*e_phi_s/sum(phi_s);
end







end