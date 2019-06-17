function [ q_stack, e_prio_stack, e_post_stack,  e_post_smooth_stack ] = ISCLP( A, Psi_w_delta, y_stack, psi_s_stack, H_stack, Psi_w_tilde_init, beta )
% [ q_stack, e_prio_stack, e_post_stack,  e_post_smooth_stack ] = ISCLP( A, Psi_w_delta, y_stack, psi_s_stack, H_stack, Psi_w_tilde_init, beta )
% runs the ISCLP Kalman filter in one frequency bin.
%
% IN:
% A                       state transition matrix - statedim x statedim
% Psi_w_delta             process noise correlation matrix - statedim x statedim
% y_stack                 microphone signal - channels x frames
% psi_s_stack             target signal PSD - frames
% H_stack                 RETF - channels (x sources) x frames (if more than one source)
% Psi_w_tilde_init        initial state estimation error correlation matrix - statedim x statedim
% beta                    gain decay limitation
%
% OUT:
% q_stack                 MF output - frames
% e_prio_stack            ISCLP output, prior - frames
% e_post_stack            ISCLP output, posterior - frames
% e_post_smooth_stack     ISCLP output, smooth posterior - frames


% get dimensions
if ismatrix(H_stack)
    N = 1;
    [M, numFrames] = size(H_stack);
else
    [M, N, numFrames] = size(H_stack);
end
stateDim = size(Psi_w_tilde_init,1);

% initKF
%
w_hat = zeros(stateDim,1);
Psi_w_tile = Psi_w_tilde_init; 
%
y_old = zeros(M,1);
%
u_LP = zeros(stateDim-M+1,1);
%
smooth_gain = 1;

% init output
q_stack = zeros(1, numFrames);
e_prio_stack = zeros(1, numFrames);
e_post_stack = zeros(1, numFrames);
e_post_smooth_stack = zeros(1, numFrames);


for i_frame = 1:numFrames
    
    %%% Load Data
    %
    y = y_stack(:,i_frame);
    psi_s = psi_s_stack(i_frame);
    if N == 1
        H = H_stack(:,i_frame);
    else
        H = H_stack(:,:,i_frame);
    end
    
    
    %%% Spatio-Temporal Pre-Processing
    %
    % compute g, B
    g = sum(H/(H'*H),2);
    Btmp = eye(M)-(H/(H'*H))*H';
    B = Btmp(:,1:M-N);
    % update q, u
    q   = g'*y;
    u_SC = B'*y;
    u_LP = [y_old; u_LP(1:end-M)];
    u = [u_SC; u_LP];
    
   
    %%% Kalman Filter
    %
    % time update: state estimate
    w_hat =  A*w_hat;
    % time update: state estimation error correlation matrix
    Psi_w_tile = A'*Psi_w_tile*A + Psi_w_delta;
    % symmetrize (just in case, avoiding accuracy issues)
    Psi_w_tile = 0.5*(Psi_w_tile+Psi_w_tile');
    % error
    e_prio_conj  = conj(q) - u'*w_hat;
    % error PSD
    psi_e = real(u'*Psi_w_tile*u) + psi_s + eps;
    % Kalman gain
    k = Psi_w_tile*u/psi_e;
    % measurement update: state estimation error correlation matrix
    w_hat = w_hat + k*e_prio_conj;
    % measurement update: state estimation error correlation matrix
    Psi_w_tile = Psi_w_tile - k*(u'*Psi_w_tile);

    
    %%% Spectral Post Processing
    %
    % compute gains
    gain        = psi_s/psi_e;
    smooth_gain = max(gain, beta*smooth_gain);
    %
    % apply gains
    e_post_conj         = gain*e_prio_conj;
    e_post_smooth_conj  = smooth_gain*e_prio_conj;


    %%% Save Data
    %
    % save
    q_stack(:,i_frame)              = q;
    e_prio_stack(:,i_frame)         = conj(e_prio_conj);
    e_post_stack(:,i_frame)         = conj(e_post_conj);
    e_post_smooth_stack(:,i_frame)  = conj(e_post_smooth_conj);
    % save previous mic. signal
    y_old = y;
    
end

end
