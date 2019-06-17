function [X, f] = calc_STFT(x, fs, win, N_STFT, R_STFT, sides)
% [X, f] = calc_STFT(x, fs, win, N_STFT, R_STFT, sides)
% performs the STFT.
%
% IN:
% x         signal - samples x channels
% fs        sampling frequency
% win       window function
% N_STFT    frame length
% R_STFT    frame shift
% sides     {'onesided', 'twosided'}, return either onesided or twosided STFT
%
% OUT:
% X         STFT tensor - freqbins x frames x channels
% f         frequency vector

% use only half FFT spectrum
N_STFT_half = N_STFT/2 + 1;

% get frequency vector
f = linspace(0,fs/2,N_STFT_half);
if strcmp(sides, 'twosided')
    f = [f, -f(end-1:-1:2)];
end

% init
L = floor((length(x) - N_STFT + R_STFT)/R_STFT);
M = size(x,2);
switch sides
    case 'onesided'
        X = zeros(N_STFT_half, L, size(x,2));   
    case 'twosided'
        X = zeros(N_STFT, L, M);   
end

for m = 1:M
    for l = 1:L % Frame index
        x_frame = x((l-1)*R_STFT+1:(l-1)*R_STFT+N_STFT, m);
        X_frame = fft(win.*x_frame);
        switch sides
            case 'onesided'
                X(:,l,m) = X_frame(1:N_STFT_half, :);
            case 'twosided'              
                X(:,l,m) = X_frame;
        end
    end
end

end