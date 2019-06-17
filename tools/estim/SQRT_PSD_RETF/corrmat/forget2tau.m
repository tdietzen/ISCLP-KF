function [ tau ] = forget2tau( forget, R_STFT, fs )
% [ tau ] = forget2tau( forget, R_STFT, fs )
% converts forgeting factor to time constant.
% 
% IN:
% forget    forgetting factor
% R_STFT    frame shift
% fs        sampling frequency
% 
% OUT:
% tau       time constant

tau = -R_STFT./(fs*log(forget));

end

