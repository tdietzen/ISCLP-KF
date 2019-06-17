function Gamma = calc_diffcoherence(micPos,N_STFT,fs,c,reg,type)
% Gamma = calc_diffcoherence(micPos,N_STFT,fs,c,reg,type)
% calculates diffuse coherence matrix.
%
% IN:
% micPos    microphone positions - channels x coordinates
% N_STFT    STFT frame length 
% fs        sampling rate
% c         speed of sound 
% reg       regularization (avoids ill-conditioned matrices at very low frequencies)
% type      {'spherical', 'cylindrical'}, coherence type
%
% OUT:
% Gamma     diffuse coherence matrix - freqbins x 1 x channels x channels

if nargin<6
    type = 'spherical';
end

N_STFT_half = N_STFT/2 + 1;
f = linspace(0,fs/2,N_STFT_half); 

M = size(micPos,1);
Gamma = (1+reg)*ones(N_STFT_half,1,M,M);

for m_out = 1:M-1
    for m_in = m_out+1:M
        d = norm(micPos(m_out,:)-micPos(m_in,:));
            switch type
                case 'spherical'
                    Gamma(:,1,m_out,m_in) = sinc(2*f*d/c);
                case 'cylindrical'
                    Gamma(:,1,m_out,m_in) =  besselj(0, 2*pi*f*d/c);
            end
        Gamma(:,1,m_in,m_out) = Gamma(:,1,m_out,m_in);
    end
end
