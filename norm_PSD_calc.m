function [f,pxx_norm] = norm_PSD(fs,signal,nfft,fc)
%This function calculates the normalized power spectral density of a given
%signal

[pxx, f] = pwelch(signal,[],[],nfft,fs);

idx = f<=fc;

%normalize PSD
df = f(2)-f(1); %Frequency resolution
f = f(idx);
pxx_norm = pxx(idx) / sum(pxx(idx)*df);

end