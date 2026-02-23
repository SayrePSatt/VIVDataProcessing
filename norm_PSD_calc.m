function [f,pxx,pxx_norm] = norm_PSD_calc(fs,signal,nfft,fc)
%This function calculates the normalized power spectral density of a given
%signal

[pxx, f] = pwelch(signal,[],[],nfft,fs);

idx = f<fc;

% size(idx)
%normalize PSD
df = f(2)-f(1); %Frequency resolution
f = f(idx);

% pxx_norm = pxx / sum(pxx*df);
% pxx_norm = pxx_norm(idx);
pxx_db = log10(pxx);

% Normalize to maximum
pxx_norm = pxx_db - max(pxx_db);
pxx_norm = pxx_norm(idx);
pxx = pxx(idx);
% pxx_norm = pxx(idx);

end