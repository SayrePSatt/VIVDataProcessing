function [f,pxx_norm] = norm_PSD_calc(fs,signal,nfft,fc)
%This function calculates the normalized power spectral density of a given
%signal

[pxx, f] = pwelch(signal,150000,100000,nfft,fs);

idx = f<=fc;

%normalize PSD
df = f(2)-f(1); %Frequency resolution
f = f(idx);
pxx_norm = pxx / sum(pxx*df);
pxx_norm = pxx_norm(idx);
% pxx_norm = pxx(idx);

end