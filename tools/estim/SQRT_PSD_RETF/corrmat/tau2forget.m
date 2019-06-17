function [ forget ] = tau2forget( tau, R_STFT, fs )
% [ forget ] = tau2forget( tau, R_STFT, fs )
% converts forgeting factor to time constant.
% 
% IN:
% tau       time constant
% R_STFT    frame shift
% fs        sampling frequency
% 
% OUT:
% forget    forgetting factor

forget = exp(-R_STFT./(fs*tau));

end

