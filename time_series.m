clear all
close all
clc

d_sph = 0.889;
[file, location] = uigetfile('*.csv');
filename = strcat(location,file);
zerofile = extractBefore(file,25);
zerofile = [zerofile '00.00.csv'];
zerofilename = strcat(location,zerofile);

cycles = 20;
f_s = 1000;
dt = 1/f_s;

data = table2array(readtable(filename));
zero = table2array(readtable(zerofilename));
time = data(:,1);
encoder = data(:,2);
encoderoffset = mean(zero(:,2));

encoder = encoder-encoderoffset;

%% Filtering
% subplot(3,1,1)
n = length(time);
fhat = fft(encoder, n); % Compute the fast fourier transform
PSD = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)

freq = 1/(dt*n)*(0:n); % Create x-qxis of frequencies in Hz
L = 1:floor(n/2); %only plot the first half of freqs
f_c = 2;   %Cutoff frequency
[b,a] = butter(4,f_c/(f_s/2),"low");

encoder_filt = filtfilt(b,a,encoder);

%% Freq Investigation
[mx,phase,f_windowed] = psdd3_sayre(f_s,encoder_filt,12000,6000,2);
f = f_windowed;
meanpwr = mean(mx,2);

nfft = 500000;
[PSD_freq, PSD_norm] = norm_PSD_calc(f_s,encoder_filt,nfft,2);

[pwr_max peak_idx] = max(meanpwr); %Finds the max power and location after taking average)

p_f = polyfit(f(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2);

f_peak = -p_f(2)/(2*p_f(1));
T = 1/f_peak;

t = time(10*T-time(1)<time<30*T(1));