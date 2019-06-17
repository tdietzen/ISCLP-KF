function H = doa2steervec(micPos, sourceAng, N_FT_half, fs, c)
% H = doa2steervec(micPos, sourceAng, N_FT_half, fs, c)
% converts DoAs in steering vectors.
%
% IN:
% micPos         microphone positions - channels x coordinates
% sourceAng      DoA angles of sources
% fs             sampling frequency
% c              speed of sound
%
% OUT:
% H              steering vectors - freqbins x 1 x channels x sources


M = size(micPos, 1);
N_src = length(sourceAng);

H = zeros(N_FT_half, 1, M, N_src);
f = (0:N_FT_half-1)*fs/(2*(N_FT_half-1));

d = sqrt(sum((micPos-repmat(micPos(1,:), [M, 1])).^2, 2));


for n = 1:N_src
    for k = 1:N_FT_half
        delay = sin(deg2rad(sourceAng(n)))*d/c;
     %   h(k,1,:,n) = exp(1i*2*pi*f(k)*delay);
        H(k,1,:,n) = exp(-1i*2*pi*f(k)*delay);
    end
end