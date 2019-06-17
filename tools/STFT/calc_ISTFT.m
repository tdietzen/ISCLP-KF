function x = calc_ISTFT(X, win, N_STFT, R_STFT, sides)
% x = calc_ISTFT(X, win, N_STFT, R_STFT, sides)
% performs the inverse STFT.
%
% IN:
% X         STFT tensor - freqbins x frames x channels
% win       window function
% N_STFT    frame length
% R_STFT    frame shift
% sides     {'onesided', 'twosided'}, return either onesided or twosided STFT
%
% OUT:
% x         signal - samples x channels

[~, L, M] = size(X);
if strcmp(sides, 'onesided')
    X = [X; conj(X(end-1:-1:2,:,:))];
end
x_frames = ifft(X, [], 1, 'symmetric');

% apply synthesis window
win = repmat(win, [1, L, M]);
x_frames = x_frames.*win;
x_frames = x_frames(1:N_STFT,:,:);

% init output
x = zeros(R_STFT*(L-1)+N_STFT, M);

% OLA processing
for l = 1:L
    sampIdx = (l-1)*R_STFT+1:(l-1)*R_STFT+N_STFT;
    x(sampIdx,:) = x(sampIdx,:) + squeeze(x_frames(:,l,:));
end