function [ q_stack, e_prio_stack, e_post_stack,  e_post_smooth_stack ] = ISCLP( A, Psi_w_delta, y_stack, psi_s_stack, H_stack, Psi_w_tilde_init, beta, complexity )
% [ q_stack, e_prio_stack, e_post_stack,  e_post_smooth_stack ] = ISCLP( A, Psi_w_delta, y_stack, psi_s_stack, H_stack, Psi_w_tilde_init, beta, complexity )
% runs the ISCLP Kalman filter in one frequency bin.
%
% IN:
% A                       state transition matrix - statedim x statedim - or scalar
% Psi_w_delta             process noise correlation matrix - statedim x statedim
% y_stack                 microphone signal - channels x frames
% psi_s_stack             target signal PSD - frames
% H_stack                 RETF - channels (x sources) x frames (if more than one source)
% Psi_w_tilde_init        initial state estimation error correlation matrix - statedim x statedim
% beta                    gain decay limitation
% complexity              {'L2M2', 'LM2', 'L2M', 'LM'}, defines complexity (default is 'L2M2')
%
% OUT:
% q_stack                 MF output - frames
% e_prio_stack            ISCLP output, prior - frames
% e_post_stack            ISCLP output, posterior - frames
% e_post_smooth_stack     ISCLP output, smooth posterior - frames

% default complexity
if nargin < 8 
    complexity = 'L2M2';
end
% get dimensions
if ismatrix(H_stack)
    N = 1;
    [M, numFrames] = size(H_stack);
else
    [M, N, numFrames] = size(H_stack);
end
stateDim = size(Psi_w_tilde_init,1);
L = ceil(stateDim/M);
N_T = L*M-stateDim;

% initKF
y_old = zeros(M,1);
u_LP = zeros(stateDim-M+1,1);
%
smooth_gain = 1;
%
switch complexity
    case 'L2M2'
        w_hat = zeros(stateDim,1);
        Psi_w_tile = Psi_w_tilde_init;
        
    case 'LM2'
        %
        % stora A and Psi_w_delta
        A_tmp = A;
        Psi_w_delta_tmp = Psi_w_delta;
        %
        % init cell arrays for L Kalman filters
        w_hat = cell(1,L);
        A = cell(1,L);
        Psi_w_delta = cell(1,L);
        Psi_w_tile = cell(1,L);
        %
        % get set of indices
        i = 1:stateDim;
        Theta = false(stateDim,L);
        for l = 1:L
            Theta(:,l) = ceil((i + N_T)/M) == l;
        end
        %
        % init Kalman filters
        for l = 1:L
            w_hat{l}        = zeros(sum(Theta(:,l)),1);
            if isscalar(A_tmp)
                A{l}        = A_tmp;
            else
                A{l}        = A_tmp(Theta(:,l), Theta(:,l));
            end
            Psi_w_delta{l}  = Psi_w_delta_tmp(Theta(:,l), Theta(:,l));
            Psi_w_tile{l}   = Psi_w_tilde_init(Theta(:,l), Theta(:,l));
        end
        
    case 'L2M'
        %
        % stora A and Psi_w_delta
        A_tmp = A;
        Psi_w_delta_tmp = Psi_w_delta;
        %
        % init cell arrays for M Kalman filters
        w_hat = cell(1,M);
        A = cell(1,M);
        Psi_w_delta = cell(1,M);
        Psi_w_tile = cell(1,M);
        %
        % get set of indices
        i = 1:stateDim;
        Xi  = false(stateDim,M);
        for m = 1:M
            Xi(:,m) = M*floor((i + N_T - m)/M) == i + N_T - m;
        end
        %
        % init Kalman filters
        for m = 1:M
            w_hat{m}        = zeros(sum(Xi(:,m)),1);
            if isscalar(A_tmp)
                A{m}        = A_tmp;
            else
                A{m}        = A_tmp(Xi(:,m), Xi(:,m));
            end
            Psi_w_delta{m}  = Psi_w_delta_tmp(Xi(:,m), Xi(:,m));
            Psi_w_tile{m}   = Psi_w_tilde_init(Xi(:,m), Xi(:,m));
        end
        
    case 'LM'
        %
        % use diagonal formulation
        w_hat       = zeros(stateDim,1);
        a           = diag(A);
        psi_w_delta = diag(Psi_w_delta);
        psi_w_tilde = diag(Psi_w_tilde_init);
end
    
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
    switch complexity
        case 'L2M2'
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
            %
        case 'LM2'
            %
            % init u
            u_tmp = u;
            u = cell(1,L);
            % init error and error PSD
            e_prio_conj = conj(q);
            psi_e = psi_s + eps;
            %
            for l = 1:L
                % get u{l}
                u{l} = u_tmp(Theta(:,l));
                %
                % time update: state estimate
                w_hat{l} =  A{l}*w_hat{l};
                % time update: state estimation error correlation matrix
                Psi_w_tile{l} = A{l}'*Psi_w_tile{l}*A{l} + Psi_w_delta{l};
                % symmetrize (just in case, avoiding accuracy issues)
                Psi_w_tile{l} = 0.5*(Psi_w_tile{l}+Psi_w_tile{l}');
                % error
                e_prio_conj  = e_prio_conj - u{l}'*w_hat{l};
                % error PSD
                psi_e = psi_e + real(u{l}'*Psi_w_tile{l}*u{l});
            end
            for l = 1:L
                % Kalman gain
                k = Psi_w_tile{l}*u{l}/psi_e;
                % measurement update: state estimation error correlation matrix
                w_hat{l} = w_hat{l} + k*e_prio_conj;
                % measurement update: state estimation error correlation matrix
                Psi_w_tile{l} = Psi_w_tile{l} - k*(u{l}'*Psi_w_tile{l});
            end
            %
        case 'L2M'
            %
            % init u
            u_tmp = u;
            u = cell(1,M);
            % init error and error PSD
            e_prio_conj = conj(q);
            psi_e = psi_s + eps;
            %
            for m = 1:M
                % get u{m}
                u{m} = u_tmp(Xi(:,m));
                %
                % time update: state estimate
                w_hat{m} =  A{m}*w_hat{m};
                % time update: state estimation error correlation matrix
                Psi_w_tile{m} = A{m}'*Psi_w_tile{m}*A{m} + Psi_w_delta{m};
                % symmetrize (just in case, avoiding accuracy issues)
                Psi_w_tile{m} = 0.5*(Psi_w_tile{m}+Psi_w_tile{m}');
                % error
                e_prio_conj  = e_prio_conj - u{m}'*w_hat{m};
                % error PSD
                psi_e = psi_e + real(u{m}'*Psi_w_tile{m}*u{m});
            end
            for m = 1:M
                % Kalman gain
                k = Psi_w_tile{m}*u{m}/psi_e;
                % measurement update: state estimation error correlation matrix
                w_hat{m} = w_hat{m} + k*e_prio_conj;
                % measurement update: state estimation error correlation matrix
                Psi_w_tile{m} = Psi_w_tile{m} - k*(u{m}'*Psi_w_tile{m});
            end
            %
        case 'LM'
            %
            % time update: state estimate
            w_hat =  a.*w_hat;
            % time update: state estimation error correlation matrix
            psi_w_tilde = a'.*psi_w_tilde.*a + psi_w_delta;
            % error
            e_prio_conj  = conj(q) - u'*w_hat;
            % error PSD
            psi_e = real(u'*(psi_w_tilde.*u)) + psi_s + eps;
            % Kalman gain
            k = psi_w_tilde.*u/psi_e;
            % measurement update: state estimation error correlation matrix
            w_hat = w_hat + k*e_prio_conj;
            % measurement update: state estimation error correlation matrix
            psi_w_tilde = psi_w_tilde - real(k.*conj(u).*psi_w_tilde);
            %
    end
        
    
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
