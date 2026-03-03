clear all
close all
clc

rho = 998;
d_sph = 0.0889;  %Diameter of Sphere
m_1k = 2.458347-(2/3)*0.0029;
m_6k = 2.458347; 
m_d = (4/3)*pi*(d_sph/2)^3*rho+rho*0.005^2*pi*d_sph/4; %Displaced mass
f_s = 1000;     %Sampling Frequency
C_A = 0.5;
bgColor = [255 255 255]/255;
load("pumpFit_freq2velo.mat");

datafolder = "E:\EFDL\viv_newstructure\";
topfolder = datafolder+"aftertare_newstructure\";

all_files = dir(topfolder);


for ii = 3:length(all_files)
    temp_config = all_files(ii).name;
    configs(ii-2) = convertCharsToStrings(temp_config(18:28));
    % distances =
end

uniq_configs = unique(configs);
uniq_configs = uniq_configs(contains(uniq_configs,test_diaratios) & contains(uniq_configs,test_distratios) & contains(uniq_configs,test_spring)); %Selects only the configurations selected for testing
uniq_configs = flip(uniq_configs);
uniq_configs = circshift(uniq_configs,1);

matching_tests = {};

plotting_color(1,:) = [0 0 0];
plotting_color(2:length(uniq_configs),:) = lines(length(uniq_configs)-1);
% plotting_color(length(uniq_configs),:) = [1 0 0];
plotting_color(2:end,:) = flipud(plotting_color(2:end,:));

marker_style = ["o"; "square"; "diamond"; "^"; "v"; ">"; "<"; "pentagram"; "hexagram";"*"];
marker_style = marker_style(1:length(uniq_configs));
marker_style = flipud(marker_style(1:length(uniq_configs)));
marker_style = circshift(marker_style,1,1);
 

for ii = 1:length(uniq_configs)
    red_velo_est = [];
    filematch = [];
    red_velo_est_temp= [];
    uniq_dist(ii) = extractBetween(uniq_configs(ii),1,3); %Extracting distance ratios
    uniq_dia(ii) = extractBetween(uniq_configs(ii),6,7); %Extracting diameter ratios
    kk = 1;
    for jj = 3:length(all_files)
        filename = all_files(jj).name;
        if contains(filename,uniq_configs(ii)) && endsWith(filename,'.dat')
            red_velo_est_temp = str2double(cell2mat(extractBetween(filename,39,42)));
            if any(red_velo_est_temp == test_redvelo)
                red_velo_est = [red_velo_est red_velo_est_temp];
                filematch = [filematch jj];
            end
        end
    end
    uniq_red_velo = unique(red_velo_est);
    for kk = 1:length(uniq_red_velo)
        uniq_red_velo(kk);
        temp_idx = find(uniq_red_velo(kk)==red_velo_est);
        matching_tests{ii,1}{kk} = filematch(temp_idx); %matching test indexing is: {configuration, spring, reduced velocity, matching indicies for each est. red. velo}
    end
end

%% Setting up files to read
test_distratios = ["000" "015" "040"];% "020" "030"];
test_diaratios = ["00" "10"];
test_spring = ["1k"];
% test_date = ["02_"];
test_redvelo = [8.5 14.5 20.5];

all_files = dir(topfolder);

for ii = 3:length(all_files)
    temp_config = all_files(ii).name;
    configs(ii-2) = convertCharsToStrings(temp_config(1:14));
    % distances =
end

uniq_configs = unique(configs);
uniq_configs = uniq_configs(contains(uniq_configs,test_diaratios) & contains(uniq_configs,test_distratios) & contains(uniq_configs,test_nums) & contains(uniq_configs,test_spring)); %Selects only the configurations selected for testing
uniq_configs = flip(uniq_configs);
uniq_configs = circshift(uniq_configs,1);
% uniq_configs = strrep(uniq_configs,test_nums,"");

for ii = 1:length(uniq_configs)
    uniq_dist = extractBetween(uniq_configs(ii),4,6); %Extracting distance ratios
    uniq_dia = extractBetween(uniq_configs(ii),9,10); %Extracting diameter ratios
    kk = 1;
    for jj = 3:length(all_files)
        filename = all_files(jj).name;
        if contains(filename,uniq_configs(ii)) && endsWith(filename,'.csv') && (contains(filename,'6.88') || contains(filename,'7.94'))
            run = convertCharsToStrings(filename);
            k_temp = str2double(extractBetween(run,13,13)); %Extracting the spring constant
            f_pump = str2double(extractBetween(run,25,29)); %Extracting Pump Speed
            if f_pump == 0
                continue
            else
                [U U_68_temp] = predict(mdl,f_pump,Alpha=0.05);%/(1.117645);
                U_68 = U-U_68_temp(1);
            end

            if k_temp == 1
                f_w = f_w_1k(1,1);
                c = c_1k;
                m = m_1k;
                k = k_1k;
                m_a = m_a_1k;
            else
                f_w = f_w_6k(1,1);
                c = c_6k;
                m = m_6k;
                k = k_6k;
                m_a = m_a_6k;
            end
            U_star = round(U/(f_w*d_sph),1); %Extracting reduced velocity
            
            data = table2array(readtable(topfolder+run));

            % zerofile = extractBefore(run,25);
            % zerofile = strcat(zerofile, '00.00.csv');
            % zero = table2array(readtable(topfolder+zerofile));
            time = data(:,1);
            encoder = data(:,2);
            encoderoffset = 0;%mean(zero(:,2));
            
            encoder = encoder-encoderoffset;

%% Filtering
cycles = 50;
f_s = 1000;
dt = 1/f_s;
% subplot(3,1,1)
n = length(time);
fhat = fft(encoder, n); % Compute the fast fourier transform
PSD = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)
freq = (1/(dt*n)*(0:n))/f_w; % Create x-qxis of frequencies in Hz
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

force_norm = 0.5*rho*(U^2)*pi*d_sph^2/4;
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
PSD_freq = PSD_freq/f_w;
[pwr_max peak_idx] = max(meanpwr); %Finds the max power and location after taking average)

p_f = polyfit(timeSeries(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2);

f_peak = -p_f(2)/(2*p_f(1));
T = 1/f_peak;

idx = (10*T-time(1))<time & time<((10+cycles)*T(1));
t = time(idx);
t = t-t(1);
encoder_filt_cont = encoder_filt(idx)/d_sph;

n = length(time);
fhat = fft(encoder, n); % Compute the fast fourier transform
PSD = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)
L = 1:floor(n/2);
freq = 1/(dt*n)*(0:n);

%PSD Investigation
[~, PSD_norm_C_y] = norm_PSD_calc(f_s,C_y,nfft,2);
[PSD_freq_filt, PSD_norm_C_vortex] = norm_PSD_calc(f_s,C_vortex,nfft,2);
PSD_freq_filt = PSD_freq_filt/f_w;

n = length(time);
fhat = fft(C_y, n); % Compute the fast fourier transform
PSD_C_y_filt = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)
fhat = fft(C_vortex, n); % Compute the fast fourier transform
PSD_C_vortex_filt = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)
L_filt = 1:floor(n/2);
freq_filt = 1/(dt*n)*(0:n);

% df = (dt*n);
% PSD_C_y_filt = PSD_C_y_filt / sum(PSD_C_y_filt*df);
% PSD_C_votex_filt = PSD_C_vortex_filt / sum(PSD_C_vortex_filt*df);
%% Plotting time history
close all
timeSeries = figure;
subplot(3,1,1)
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
% timeSeries.Position = [100 100 500 250];
plot(t/T,encoder_filt_cont,'k-')
xlim([0 cycles])
limits = ceil(max(abs(encoder_filt_cont))/0.5)*0.5;
% limits = 0.1;
ylim([-limits limits])
yticks([-limits 0 limits])
% xlabel('$t/T$')
ylabel('$y/D$')
distance = str2double(char(uniq_dist))/10;
if distance == 0
    L_star = 'Isolated ';
else
    L_star = ['$L^*=$' +num2str(distance) ' '];
end
title([L_star ' ' num2str(k_temp) 'k Config ' '$U^*=$' num2str(U_star)])
set(gca,'xticklabel',[])
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);

% figurename = strcat(extractBefore(run,12),'_',string(k_temp),'k_Ustar_', num2str(U_star*10), '_timeseries');
% exportgraphics(timeSeries,strcat('figures\timeSeries\', figurename, '.png'),'Resolution',300,'BackgroundColor',bgColor);
% saveas(timeSeries,strcat('figures\timeSeries\', figurename, '.eps'),'epsc');

%Phase time series plot
subplot(3,1,2)
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
% timeSeries.Position = [100 100 500 250];
plot(t/T,totalphase(idx),'k-')
yline(mean(totalphase),'r--')
yline(90,'b:')
xlim([0 cycles])
limits = 180;%ceil(max(abs(encoder_filt_cont))/0.5)*0.5;
% limits = 0.1;
ylim([-90 270])
yticks([-90 0 90 180 270])
% xlabel('$t/T$')
ylabel('$\phi_{total}$')
% distance = str2double(char(uniq_dist))/10;
% if distance == 0
%     L_star = 'Isolated ';
% else
%     L_star = ['$L^*=$' +num2str(distance) ' '];
% end
% title([L_star '$U^*=$' num2str(U_star)])
set(gca,'xticklabel',[])
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);

%Vortex Phase Time Series
subplot(3,1,3)
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
timeSeries.Position = [100 100 500 400];
plot(t/T,vortexphase(idx),'k-')
yline(mean(vortexphase),'r--')
yline(90,'b:')
xlim([0 cycles])
limits = 180;%ceil(max(abs(encoder_filt_cont))/0.5)*0.5;
% limits = 0.1;
ylim([-90 270])
yticks([-90 0 90 180 270])
xlabel('$t/T$')
ylabel('$\phi_{vortex}$')
% distance = str2double(char(uniq_dist))/10;
% if distance == 0
%     L_star = 'Isolated ';
% else
%     L_star = ['$L^*=$' +num2str(distance) ' '];
% end
% title([L_star '$U^*=$' num2str(U_star)])
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);

figurename = strcat(extractBefore(run,12),'_',string(k_temp),'k_Ustar_', num2str(U_star*10), '_timeseries');
exportgraphics(timeSeries,strcat('figures\timeSeries\', figurename, '.png'),'Resolution',300,'BackgroundColor',bgColor);
saveas(timeSeries,strcat('figures\timeSeries\', figurename, '.eps'),'epsc');

% % Plotting Lissajous Curves
% lissajous_fig = figure;
% 
% x = encoder_filt(1000:end-1000)/d_sph;
% y = C_y(1000:end-1000);
% 
% plot(x,y,'LineWidth',1.5,'Color','k')
% xlabel('$y/D$')
% ylabel('$C_y$')
% title([L_star '$U^*=$' num2str(U_star)])
% xbounds = max(x+0.1);
% ybounds = max(y+0.1);
% xlim([-xbounds xbounds]);
% ylim([-ybounds ybounds]);
% axis square
% set(gcf, 'color', bgColor);
% set(gca, 'color', bgColor);
% figsize = get(gcf,'Position');
% figurename = strcat(extractBefore(run,12),'_',string(k_temp),'k_Ustar_', num2str(U_star*10), '_lissajous');
% exportgraphics(lissajous_fig,strcat('figures\timeSeries\', figurename, '.png'),'Resolution',300,'BackgroundColor',bgColor);
% saveas(lissajous_fig,strcat('figures\timeSeries\', figurename, '.eps'),'epsc');

%% Frequency and PSD Plotting
PSD_figure = figure;
tl = tiledlayout(3,1);
nexttile
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
% timeSeries.Position = [100 100 500 250];
semilogy(freq(L)/f_w,PSD(L),'k-')
hold on
semilogy(PSD_freq,PSD_norm,'r-')
% xline(f_peak,'r-')
% xline(St*U/d_sph,'b:')
xlim([0 5])
% limits = ceil(max(abs(encoder_filt_cont))/0.5)*0.5;
% % limits = 0.1;
% ylim([-limits limits])
% yticks([-limits 0 limits])
% xlabel('$t/T$')
ylabel('$x_{PSD}$')
distance = str2double(char(uniq_dist))/10;
if distance == 0
    L_star = 'Isolated ';
else
    L_star = ['$L^*=$' +num2str(distance) ' '];
end
title([L_star ' ' num2str(k_temp) 'k Config ' '$U^*=$' num2str(U_star)])
set(gca,'xticklabel',[])
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);

% figurename = strcat(extractBefore(run,12),'_',string(k_temp),'k_Ustar_', num2str(U_star*10), '_timeseries');
% exportgraphics(timeSeries,strcat('figures\timeSeries\', figurename, '.png'),'Resolution',300,'BackgroundColor',bgColor);
% saveas(timeSeries,strcat('figures\timeSeries\', figurename, '.eps'),'epsc');

%Phase time series plot
nexttile
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
% timeSeries.Position = [100 100 500 250];
semilogy(freq_filt(L)/f_w,PSD_C_y_filt(L),'k-')
hold on
semilogy(PSD_freq_filt,PSD_norm_C_y,'r-')
% xline(f_peak,'r-')
% xline(St*U/d_sph,'b:')
xlim([0 5])
% limits = ceil(max(abs(encoder_filt_cont))/0.5)*0.5;
% % limits = 0.1;
% ylim([-limits limits])
% yticks([-limits 0 limits])
% xlabel('$t/T$')
ylabel('$C_{y,PSD}$')
distance = str2double(char(uniq_dist))/10;
if distance == 0
    L_star = 'Isolated ';
else
    L_star = ['$L^*=$' +num2str(distance) ' '];
end
% title([L_star '$U^*=$' num2str(U_star)])
set(gca,'xticklabel',[])
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);

%Vortex Phase Time Series
nexttile
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
% timeSeries.Position = [100 100 500 250];
semilogy(freq_filt(L)/f_w,PSD_C_vortex_filt(L),'k-')
hold on
semilogy(PSD_freq_filt,PSD_norm_C_vortex,'r-')
% xline(f_peak,'r-')
% xline(St*U/d_sph,'b:')
xlim([0 5])
% limits = ceil(max(abs(encoder_filt_cont))/0.5)*0.5;
% % limits = 0.1;
% ylim([-limits limits])
% yticks([-limits 0 limits])
xlabel('$f^*$')
ylabel('$C_{vortex,PSD}$')
distance = str2double(char(uniq_dist))/10;
if distance == 0
    L_star = 'Isolated ';
else
    L_star = ['$L^*=$' +num2str(distance) ' '];
end
% title([L_star '$U^*=$' num2str(U_star)])
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);

position_temp = tl.InnerPosition;
norm_bottom = position_temp(2);
norm_top = position_temp(2)+position_temp(4);
ax_offset = get(tl,'Position');
ax = gca;
ax_limits = get(ax,'XLim');
axis_upper = ax_limits(2);
axis_lower = ax_limits(1);
St_norm_value = ax_offset(1)+((St*U/d_sph/f_w)-axis_lower)/(axis_upper-axis_lower)*ax_offset(3);
fstar_norm_value = ax_offset(1)+((f_peak/f_w)-axis_lower)/(axis_upper-axis_lower)*ax_offset(3);

annotation(PSD_figure, 'line', [St_norm_value St_norm_value], [norm_bottom norm_top], 'Color', 'b', 'LineWidth', 1.5,'LineStyle',':');
annotation(PSD_figure, 'line', [fstar_norm_value fstar_norm_value], [norm_bottom norm_top], 'Color', 'r', 'LineWidth', 1.5,'LineStyle','--');

figurename = strcat(extractBefore(run,12),'_',string(k_temp),'k_Ustar_', num2str(U_star*10), '_spectral');
exportgraphics(PSD_figure,strcat('figures\spectralAnalysis\', figurename, '.png'),'Resolution',300,'BackgroundColor',bgColor);
saveas(PSD_figure,strcat('figures\spectralAnalysis\', figurename, '.eps'),'epsc');
        end
    end
end

% [file, location] = uigetfile('*.csv');
% filename = strcat(location,file);
% zerofile = extractBefore(file,25);
% zerofile = [zerofile '00.00.csv'];
% zerofilename = strcat(location,zerofile);
% load("pumpFit_freq2velo.mat");
% decay_location = strcat(extractBefore(location,'testData\'),'freeDecay\');
% %%
% 
% 
% 
% pump_f = str2double(extractBetween(file,25,29)); %extracting pump frequency
% pump_U = predict(mdl,pump_f);
% 
% % test_num = extractBefore(file,3);
% spring_num = char(extractBetween(file,13,13));
% distance = str2double(char(extractBetween(file,4,6)))/10;
% 
% freedecay_location = strcat(decay_location,'latestFreedecay\freedecay_',spring_num,'k_water.dat');
% f_w = table2array(readtable(freedecay_location));
% f_w = f_w(1,1);
% freedecay_location = strcat(decay_location,'latestFreedecay\freedecay_',spring_num,'k_air.dat');
% f_n = table2array(readtable(freedecay_location));
% zeta = f_n(1,2);
% f_n = f_n(1,1);
% 
% m_a = ((f_n/f_w).^2-1)*m;
% omegana = 2*pi*f_n;
% k = m*omegana.^2;
% c = zeta*2*sqrt(m*k);
% 
% U_star = round(pump_U/(f_w*d_sph),1);
% 
% data = table2array(readtable(filename));
% zero = table2array(readtable(zerofilename));
% time = data(:,1);
% encoder = data(:,2);
% encoderoffset = mean(zero(:,2));
% 
% encoder = encoder-encoderoffset;
% 
% %% Filtering
% cycles = 50;
% f_s = 1000;
% dt = 1/f_s;
% % subplot(3,1,1)
% n = length(time);
% fhat = fft(encoder, n); % Compute the fast fourier transform
% PSD = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)
% 
% freq = 1/(dt*n)*(0:n); % Create x-qxis of frequencies in Hz
% L = 1:floor(n/2); %only plot the first half of freqs
% f_c = 2;   %Cutoff frequency
% [b,a] = butter(4,f_c/(f_s/2),"low");
% 
% encoder_filt = encoder;%filtfilt(b,a,encoder)/d_sph;
% 
% velo = FivePointDiff(encoder_filt,f_s)';
% velo = filtfilt(b,a,velo);
% acc = FivePointDiff(velo,f_s)';
% acc = filtfilt(b,a,acc);
% 
% data_length = length(acc);
% time = time(1:data_length);
% encoder_filt = encoder_filt(1:data_length);
% velo =  velo(1:data_length);
% 
% F = m.*acc+c.*velo+k*encoder_filt;
% F_pot = -m_a*acc;%-C_A*m_d*acc;
% F_vortex = F - F_pot;
% 
% force_norm = 0.5*rho*(pump_U^2)*pi*d_sph^2/4;
% C_y = F/force_norm;
% C_pot = F_pot/force_norm;
% C_vortex = F_vortex/force_norm;
% 
% [totalphase, vortexphase] = retrievephase1(encoder_filt,C_y,C_vortex);
% %% Freq Investigation
% [mx,phase,f_windowed] = psdd3_sayre(f_s,encoder_filt,12000,6000,2);
% timeSeries = f_windowed;
% meanpwr = mean(mx,2);
% 
% nfft = 500000;
% [PSD_freq, PSD_norm] = norm_PSD_calc(f_s,encoder_filt,nfft,2);
% 
% [pwr_max peak_idx] = max(meanpwr); %Finds the max power and location after taking average)
% 
% p_f = polyfit(timeSeries(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2);
% 
% f_peak = -p_f(2)/(2*p_f(1));
% T = 1/f_peak;
% 
% idx = (10*T-time(1))<time & time<((10+cycles)*T(1));
% t = time(idx);
% t = t-t(1);
% encoder_filt_cont = encoder_filt(idx)/d_sph;
% %% Plotting time history
% plot_color = [255 255 255]/255;
% close all
% timeSeries = figure;
% set(gcf, 'color', plot_color);
% set(gca, 'color', plot_color);
% timeSeries.Position = [100 100 500 250];
% plot(t/T,encoder_filt_cont,'k-')
% xlim([0 cycles])
% limits = ceil(max(abs(encoder_filt_cont))/0.5)*0.5;
% % limits = 0.1;
% ylim([-limits limits])
% yticks([-limits 0 limits])
% xlabel('$t/T$')
% ylabel('$y/D$')
% if distance == 0
%     L_star = 'Isolated ';
% else
%     L_star = ['$L^*=$' +num2str(distance) ' '];
% end
% title([L_star '$U^*=$' num2str(U_star)])
% set(gcf, 'color', plot_color);
% set(gca, 'color', plot_color);
% 
% figurename = [extractBefore(file,11) '_' num2str(U_star*10) '_timeseries'];
% exportgraphics(timeSeries,['figures\' figurename '.png'],'Resolution',300,'BackgroundColor',[255 255 255]/255);
% saveas(timeSeries,['figures\' figurename '.eps'],'epsc');
% 
% %% Plotting Lissajous Curves
% lissajous_fig = figure;
% 
% x = encoder_filt(1000:end-1000)/d_sph;
% y = C_y(1000:end-1000);
% 
% plot(x,y,'LineWidth',3,'Color','k')
% xlabel('$y/D$')
% ylabel('$C_y$')
% title([L_star '$U^*=$' num2str(U_star)])
% xbounds = max(x+0.1);
% ybounds = max(y+0.1);
% xlim([-xbounds xbounds]);
% ylim([-ybounds ybounds]);
% axis square
% set(gcf, 'color', plot_color);
% set(gca, 'color', plot_color);
% figsize = get(gcf,'Position');
% figurename = [extractBefore(file,11) '_' num2str(U_star*10) '_lissajous'];
% exportgraphics(lissajous_fig,['figures\' figurename '.png'],'Resolution',300,'BackgroundColor',[255 255 255]/255);
% saveas(lissajous_fig,['figures\' figurename '.eps'],'epsc');

%% Animated plots

% lissajous_ani_fig = figure;
% lissajous_ani_fig.Position = figsize;
% axis square
% box on
% xlim([-xbounds xbounds]);
% ylim([-ybounds ybounds]);
% title(['$U^*=$' num2str(U_star) ' $L^*=$' num2str(distance)])
% xlabel('$y/D$')
% ylabel('$C_y$')
% set(gcf, 'color', plot_color);
% set(gca, 'color', plot_color);
% 
% lissajous_ani = animatedline('Color','k','LineWidth',3,'MaximumNumPoints',20000);
% lissajous_ani_point = animatedline('Color','r','LineStyle','none','MaximumNumPoints',1,'Marker','o','MarkerSize',12,'MarkerFaceColor','b');
% numpoints=100000;
% 
% filename = 'figures\lissajous.gif';
% framerate = 60;
% 
% for k = 1:50:50*700
%     xvec = x(k:k+99);
%     yvec = y(k:k+99);
%     addpoints(lissajous_ani, xvec, yvec);
%     addpoints(lissajous_ani_point, xvec(end), yvec(end));
%     drawnow
%     % Capture the plot as an image
%     frame = getframe(lissajous_ani_fig);
%     img = frame2im(frame);
%     [imind, cm] = rgb2ind(img, 256);
%     % Write to the GIF File
%     if k == 1
%         imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1/framerate); % Adjust DelayTime as needed
%     else
%         imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/framerate);
%     end
% end