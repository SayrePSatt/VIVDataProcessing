function [f,pxx_norm] = FFT_PSD_calc(fs,signal,nfft,fc)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
N = length(signal);
signal_fft = fft(signal);
f = (0:N-1)*(fs/N);

PSD = (1/(fs*N)) * abs(signal_fft).^2;
PSD = PSD(1:N/2+1);
PSD(2:end-1) = 2*PSD(2:end-1); % Only positive freqs, double except endpoints
f = f(1:N/2+1);

idx = f<fc;
f = f(idx);
% size(idx)
%normalize PSD
% df = f(2)-f(1); %Frequency resolution
% f = f(idx);

% pxx_norm = pxx / sum(pxx*df);
% pxx_norm = pxx_norm(idx);
pxx_db = log10(PSD(idx));

% Normalize to maximum
pxx_norm = pxx_db - max(pxx_db);
pxx_norm = pxx_norm(idx);
PSD = PSD(idx);
% pxx_norm = pxx(idx);


% smooth_PSD = sgolayfilt(PSD,1,5);
end