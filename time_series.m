clear all
close all
clc

rho = 998;
d_sph = 0.0889;  %Diameter of Sphere
m = 2.42947;    %Oscillating Mass. 2.4295 for 90mm setup, 1.916 for 80mm setup
m_d = (4/3)*pi*(d_sph/2)^3*rho+rho*0.005^2*pi*d_sph/4; %Displaced mass
f_s = 1000;     %Sampling Frequency
C_A = 0.5;

[file, location] = uigetfile('*.csv');
filename = strcat(location,file);
zerofile = extractBefore(file,25);
zerofile = [zerofile '00.00.csv'];
zerofilename = strcat(location,zerofile);
load("pumpFit_freq2velo.mat");
decay_location = strcat(extractBefore(location,'testData\'),'freeDecay\');
%%

pump_f = str2double(extractBetween(file,25,29)); %extracting pump frequency
pump_U = predict(mdl,pump_f);

test_num = extractBefore(file,3);
spring_num = char(extractBetween(file,13,13));
distance = str2double(char(extractBetween(file,4,6)))/10;

freedecay_location = strcat(decay_location,test_num,'_test\freedecay_',spring_num,'k_water.dat');
f_w = table2array(readtable(freedecay_location));
f_w = f_w(1,1);
freedecay_location = strcat(decay_location,test_num,'_test\freedecay_',spring_num,'k_air.dat');
f_n = table2array(readtable(freedecay_location));
zeta = f_n(1,2);
f_n = f_n(1,1);

m_a = ((f_n/f_w).^2-1)*m;
omegana = 2*pi*f_n;
k = m*omegana.^2;
c = zeta*2*sqrt(m*k);

U_star = round(pump_U/(f_w*d_sph),1);

data = table2array(readtable(filename));
zero = table2array(readtable(zerofilename));
time = data(:,1);
encoder = data(:,2);
encoderoffset = mean(zero(:,2));

encoder = encoder-encoderoffset;

%% Filtering
cycles = 50;
f_s = 1000;
dt = 1/f_s;
% subplot(3,1,1)
n = length(time);
fhat = fft(encoder, n); % Compute the fast fourier transform
PSD = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)

freq = 1/(dt*n)*(0:n); % Create x-qxis of frequencies in Hz
L = 1:floor(n/2); %only plot the first half of freqs
f_c = 2;   %Cutoff frequency
[b,a] = butter(4,f_c/(f_s/2),"low");

encoder_filt = encoder;%filtfilt(b,a,encoder)/d_sph;

velo = FivePointDiff(encoder_filt,f_s)';
velo = filtfilt(b,a,velo);
acc = FivePointDiff(velo,f_s)';
acc = filtfilt(b,a,acc);

data_length = length(acc);
time = time(1:data_length);
encoder_filt = encoder_filt(1:data_length);
velo =  velo(1:data_length);

F = m.*acc+c.*velo+k*encoder_filt;
F_pot = -m_a*acc;%-C_A*m_d*acc;
F_vortex = F - F_pot;

force_norm = 0.5*rho*(pump_U^2)*pi*d_sph^2/4;
C_y = F/force_norm;
C_pot = F_pot/force_norm;
C_vortex = F_vortex/force_norm;

[totalphase, vortexphase] = retrievephase1(encoder_filt,C_y,C_vortex);
%% Freq Investigation
[mx,phase,f_windowed] = psdd3_sayre(f_s,encoder_filt,12000,6000,2);
timeSeries = f_windowed;
meanpwr = mean(mx,2);

nfft = 500000;
[PSD_freq, PSD_norm] = norm_PSD_calc(f_s,encoder_filt,nfft,2);

[pwr_max peak_idx] = max(meanpwr); %Finds the max power and location after taking average)

p_f = polyfit(timeSeries(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2);

f_peak = -p_f(2)/(2*p_f(1));
T = 1/f_peak;

idx = (10*T-time(1))<time & time<((10+cycles)*T(1));
t = time(idx);
t = t-t(1);
encoder_filt_cont = encoder_filt(idx)/d_sph;
%% Plotting time history
plot_color = [238 238 238]/255;
close all
timeSeries = figure;
set(gcf, 'color', plot_color);
set(gca, 'color', plot_color);
timeSeries.Position = [100 100 500 250];
plot(t/T,encoder_filt_cont,'k-')
xlim([0 cycles])
limits = ceil(max(abs(encoder_filt_cont))/0.5)*0.5;
% limits = 0.1;
ylim([-limits limits])
yticks([-limits 0 limits])
xlabel('$t/T$')
ylabel('$y/D$')
title(['$U^*=$' num2str(U_star) ' $L^*=$' num2str(distance)])
set(gcf, 'color', plot_color);
set(gca, 'color', plot_color);

figurename = [extractBefore(file,11) '_' num2str(U_star*10) '_timeseries'];
exportgraphics(timeSeries,['figures\' figurename '.png'],'Resolution',300,'BackgroundColor',[238 238 238]/255);
saveas(timeSeries,['figures\' figurename '.eps'],'epsc');

%% Plotting Lissajous Curves
lissajous_fig = figure;

x = encoder_filt(1000:end-1000)/d_sph;
y = C_y(1000:end-1000);

plot(x,y,'LineWidth',3,'Color','k')
xlabel('$y/D$')
ylabel('$C_y$')
title(['$U^*=$' num2str(U_star) ' $L^*=$' num2str(distance)])
xbounds = max(x+0.1);
ybounds = max(y+0.1);
xlim([-xbounds xbounds]);
ylim([-ybounds ybounds]);
axis square
set(gcf, 'color', plot_color);
set(gca, 'color', plot_color);
figsize = get(gcf,'Position');
figurename = [extractBefore(file,11) '_' num2str(U_star*10) '_lissajous'];
exportgraphics(lissajous_fig,['figures\' figurename '.png'],'Resolution',300,'BackgroundColor',[238 238 238]/255);
saveas(lissajous_fig,['figures\' figurename '.eps'],'epsc');

%% Animated plots

lissajous_ani_fig = figure;
lissajous_ani_fig.Position = figsize;
axis square
box on
xlim([-xbounds xbounds]);
ylim([-ybounds ybounds]);
title(['$U^*=$' num2str(U_star) ' $L^*=$' num2str(distance)])
xlabel('$y/D$')
ylabel('$C_y$')
set(gcf, 'color', plot_color);
set(gca, 'color', plot_color);

lissajous_ani = animatedline('Color','k','LineWidth',3,'MaximumNumPoints',20000);
lissajous_ani_point = animatedline('Color','r','LineStyle','none','MaximumNumPoints',1,'Marker','o','MarkerSize',12,'MarkerFaceColor','b');
numpoints=100000;

filename = 'figures\lissajous.gif';
framerate = 60;

for k = 1:50:50*700
    xvec = x(k:k+99);
    yvec = y(k:k+99);
    addpoints(lissajous_ani, xvec, yvec);
    addpoints(lissajous_ani_point, xvec(end), yvec(end));
    drawnow
    % Capture the plot as an image
    frame = getframe(lissajous_ani_fig);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);
    % Write to the GIF File
    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1/framerate); % Adjust DelayTime as needed
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/framerate);
    end
end