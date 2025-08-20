clear all
close all
clc

d_sph = 0.0889;
[file, location] = uigetfile('*.csv');
filename = strcat(location,file);
zerofile = extractBefore(file,25);
zerofile = [zerofile '00.00.csv'];
zerofilename = strcat(location,zerofile);
load("pumpFit_freq2velo.mat");

pump_f = str2double(extractBetween(file,25,29)); %extracting pump frequency
pump_U = predict(mdl,pump_f);

test_num = extractBefore(file,3);
spring_num = char(extractBetween(file,13,13));

freedecay_location = strcat(location,test_num,'_test\freedecay_',spring_num,'k_water.dat');
f_w = table2array(readtable(freedecay_location));
f_w = f_w(1,1);

U_star = round(pump_U/(f_w*d_sph),1);

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

encoder_filt = filtfilt(b,a,encoder)/d_sph;

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

idx = (10*T-time(1))<time & time<((10+cycles)*T(1));
t = time(idx);
t = t-t(1);
encoder_filt = encoder_filt(idx);
%% Plotting
close all
f = figure;
f.Position = [100 100 500 250];
figure(f)
plot(t/T,encoder_filt,'k-')
xlim([0 cycles])
limits = round(max(abs(encoder_filt))/0.5)*0.5;
% limits = 0.1;
ylim([-limits limits])
yticks([-limits 0 limits])
xlabel('$t/T$')
ylabel('$y/D$')
title(['$U^*=$' num2str(U_star)])

figurename = [extractBefore(file,11) '_' num2str(U_star*10) '_timeseries'];
exportgraphics(f,['figures\' figurename '.jpg'],'Resolution',300);
saveas(f,['figures\' figurename '.eps'],'epsc');