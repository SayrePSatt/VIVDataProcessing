%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: ReducedFreq_Response.m
% Version: 1
% Date: 04/28/2025
% Author: Sayre Satterwhite (sayreps@umich.edu)
% Description: Takes displacement data for a sphere in VIV and determines
% relavant quantities, compares with historical results
% This is adapted from a previous version for use with new data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all %#ok<CLALL>
close all
clc

%% Options for plotting
plot_legends = 1; %0 to not plot legends, 1 to plot legends
plot_reference = 0; %0 to not plot references
plot_errors = 1; %0 to not plot errorbars
single_test = 1; %Use for plotting the spectrogram curves and mean peaks curve

test_distratios = ["000" "015" "040" "100"];% "020" "025" "030" "040" "050" "060" "070" "100"];% "020" "030"];
test_diaratios = ["00" "10"]; %"06" "08"];
test_spring = ["1k"];%["6k" "1k"];

bgColor = [255 255 255]/255;

%% Experiment Specification
datafolder = "D:\vivscratch_complete\";
topfolder = datafolder+"aftertare\";

rho = 998;
d_sph = 0.0889;  %Diameter of Sphere
m_1k = 2.458347-(2/3)*0.0029;
m_6k = 2.458347;    %Oscillating Mass. 2.4295 for 90mm setup, 1.916 for 80mm setup
m_d = (4/3)*pi*(d_sph/2)^3*rho+rho*0.005^2*pi*d_sph/4; %Displaced mass
f_s = 1000;     %Sampling Frequency
C_A = 0.5;     %Added mass coefficient
temp_1k = table2array(readtable(datafolder+"freeDecay/1k_09_26_2025/freedecay_1k_air.dat"));
f_n_1k(1,:) = temp_1k(1,:);
f_n_1k_95 = temp_1k(2,1);
% zeta_1k_95 = temp_1k(2,2);
temp_1k = table2array(readtable(datafolder+"freeDecay/1k_09_26_2025/freedecay_1k_water.dat"));
f_w_1k_95 = temp_1k(2,1);
zeta_1k_95 = temp_1k(2,2);
f_w_1k(1,:) = temp_1k(1,:);
f_n_1k(2,:) = f_n_1k(1,:);
f_w_1k(2,:) = f_w_1k(1,:);
f_n_1k(3,:) = f_n_1k(1,:);
f_w_1k(3,:) = f_w_1k(1,:);
f_n_1k(4,:) = f_n_1k(1,:);
f_w_1k(4,:) = f_w_1k(1,:);

% temp_1k = table2array(readtable(datafolder+"freeDecay/1k_08_18_2025/freedecay_1k_air.dat"));
% f_n_1k(5,:) = temp_1k(1,:);
% temp_1k = table2array(readtable(datafolder+"freeDecay/1k_08_18_2025/freedecay_1k_water.dat"));
% f_w_1k(5,:) = temp_1k(1,:);
% f_n_1k(6,:) = f_n_1k(5,:);
% f_w_1k(6,:) = f_w_1k(5,:);
% f_n_1k(7,:) = f_n_1k(5,:);
% f_w_1k(7,:) = f_w_1k(5,:);
% f_n_1k(8,:) = f_n_1k(5,:);
% f_w_1k(8,:) = f_w_1k(5,:);

temp_6k = table2array(readtable(datafolder+"freeDecay/6k_09_26_2025/freedecay_6k_air.dat"));
f_n_6k(1,:) = temp_6k(1,:);
f_n_6k_95 = temp_6k(2,1);
% zeta_6k_95 = temp_1k(2,2);
temp_6k = table2array(readtable(datafolder+"freeDecay/6k_09_26_2025/freedecay_6k_water.dat"));
f_w_6k_95 = temp_6k(2,1);
zeta_6k_95 = temp_6k(2,2);
f_w_6k(1,:) = temp_6k(1,:);
f_w_6k(2,:) = f_w_6k(1,:);
f_w_6k(3,:) = f_w_6k(1,:);
f_w_6k(4,:) = f_w_6k(1,:);
f_n_6k(2,:) = f_n_6k(1,:);
f_n_6k(3,:) = f_n_6k(1,:);
f_n_6k(4,:) = f_n_6k(1,:);

% temp_6k = table2array(readtable(datafolder+"freeDecay/6k_08_18_2025/freedecay_6k_air.dat"));
% f_n_6k(5,:) = temp_6k(1,:);
% temp_6k = table2array(readtable(datafolder+"freeDecay/6k_08_18_2025/freedecay_6k_water.dat"));
% f_w_6k(5,:) = temp_6k(1,:);
% f_w_6k(6,:) = f_w_6k(5,:);
% f_w_6k(7,:) = f_w_6k(5,:);
% f_w_6k(8,:) = f_w_6k(5,:);
% f_n_6k(6,:) = f_n_6k(5,:);
% f_n_6k(7,:) = f_n_6k(5,:);
% f_n_6k(8,:) = f_n_6k(5,:);

% f_w_1k(:,1) = 0.2413;
% f_n_1k(:,1) = 0.2495;
% f_1_1k(:,2) = 0.0068;
% 
% f_w_6k(:,1) = 0.5697;
% f_n_6k(:,1) = 0.5917;
% f_n_6k(:,2) = 0.0029;

m_a_1k = ((f_n_1k(:,1)./f_w_1k(:,1)).^2-1)*m_1k; %test
m_a_6k = ((f_n_6k(:,1)./f_w_6k(:,1)).^2-1)*m_6k;

% ca_1k = m_a_1k/m_d
% ca_6k = m_a_6k/m_d

% m_a_6k(:) = m_d/2;
% m_a_1k(:) = m_d/2;
St = 0.19;
St_68 = 0.005;
omegana_1k = 2*pi*f_n_1k(:,1);
k_1k = m_1k*omegana_1k.^2; %5.375; %(f_n(1)*2*pi)^2*m
% k_1k(:) = 6;
% k_1k = k_1k-0.6;
omegana_6k = 2*pi*f_n_6k(:,1);
k_6k = m_6k*omegana_6k.^2;
% k_6k = k_6k - 1.0;
% k_6k(:) = 33.7;

c_1k = 4*pi*f_n_1k(:,2).*m_1k.*f_n_1k(:,1);
c_6k = 4*pi*f_n_6k(:,2).*m_6k.*f_n_6k(:,1);
m_star_1k = m_1k/m_d;
m_star_6k = m_6k/m_d;
mass_damp_1k = (m_star_1k+C_A)*f_n_1k(1,2);
mass_damp_6k = (m_star_6k+C_A)*f_n_6k(1,2);
% scruton = 2*m*f_n_1k(2)/(rho*d_sph^2); 

load("pumpFit_freq2velo.mat");

diagnose = false;

markers = ['s' 'd' '*'];
%% Importing external data for comparison
sareen = csvread('sareen2018b_ampphase.csv',3); %#ok<CSVRD>
sareen(sareen==0) = NaN;
u_red_A_star_sareen = sareen(:,1);
A_star_sareen = sareen(:,2);
u_red_vortexphase_sareen = sareen(:,3);
vortexphase_sareen = sareen(:,4);
u_red_totalphase_sareen = sareen(:,5);
totalphase_sareen = sareen(:,6);

govwill = csvread('gov_2005_amplift.csv',3); %#ok<CSVRD>
govwill(govwill==0) = NaN;
u_red_A_star_govwill = govwill(:,9);
A_star_govwill = govwill(:,10);
u_red_totalforce_govwill = govwill(:,1);
totalforce_govwill = govwill(:,2);
u_red_totalphase_govwill = govwill(:,3);
totalphase_govwill = govwill(:,4);
u_red_vortexforce_govwill = govwill(:,5);
vortexforce_govwill = govwill(:,6);
u_red_vortexphase_govwill = govwill(:,7);
vortexphase_govwill = govwill(:,8);

unorm_reference = csvread('unorm_reference.csv',3); %#ok<CSVRD>
unorm_reference(unorm_reference==0) = NaN;
unorm_unorm_sareen = unorm_reference(:,1);
unorm_Astar_sareen = unorm_reference(:,2);
unorm_unorm_govwill = rmmissing(unorm_reference(:,3));
unorm_Astar_govwill = rmmissing(unorm_reference(:,4));

griffin_plot = csvread('griffin_govwill.csv',3); %#ok<CSVRD>
griffin_plot(griffin_plot==0) = NaN;
griffin_massdamp_samp = griffin_plot(:,1);
griffin_Astar_samp = griffin_plot(:,2);
griffin_massdamp_fit = rmmissing(griffin_plot(:,3));
griffin_Astar_fit = rmmissing(griffin_plot(:,4));
%% Setting up plots for later

A_y_star_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

A_y_norm_fig = figure;
A_y_norm_fig.Position = [100 100 600 400];
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

phase_subplot_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

total_force_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

vortex_force_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

total_phase_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

vortex_phase_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

pdicy_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

f_star_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

test_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

griffin_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

if single_test == 1
    A_y_star_pctile_fig = figure;
    set(gcf, 'color', bgColor);
    set(gca, 'color', bgColor);
    hold on
    
    freq_contour_fig = figure;
    set(gcf, 'color', bgColor);
    set(gca, 'color', bgColor);
    hold on;
end
%% Setting up folder directories
all_files = dir(topfolder);


for ii = 3:length(all_files)
    temp_config = all_files(ii).name;
    configs(ii-2) = convertCharsToStrings(temp_config(4:14));
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
    uniq_dist(ii) = extractBetween(uniq_configs(ii),1,3); %Extracting distance ratios
    uniq_dia(ii) = extractBetween(uniq_configs(ii),6,7); %Extracting diameter ratios
    kk = 1;
    for jj = 3:length(all_files)
        filename = all_files(jj).name;
        if contains(filename,uniq_configs(ii)) && endsWith(filename,'.csv')
            matching_tests{ii,1}(kk) = convertCharsToStrings(filename);
            
            matching_tests{ii,2}(kk) = str2double(extractBetween(matching_tests{ii,1}(kk),13,13)); %Extracting the spring constant
            matching_tests{ii,3}(kk) = str2double(extractBetween(matching_tests{ii,1}(kk),25,29)); %Extracting Pump Speed
            f_pump = matching_tests{ii,3}(kk);
            if f_pump == 0
                U = 0.0;
                U_68 = 0;
            else
                [U U_68_temp] = predict(mdl,f_pump,Alpha=0.05);%/(1.117645);
                U_68 = U-U_68_temp(1);
            end
            matching_tests{ii,6}(kk) = U; %Extracting flow velocity
            matching_tests{ii,4}(kk) = str2double(extractBetween(matching_tests{ii,1}(kk),1,2)); %Extracting the test number
           
            k_temp = matching_tests{ii,2}(kk);
            if k_temp == 1
                f_w = f_w_1k(matching_tests{ii,4}(kk),1)
                m = m_1k;
            else
                f_w = f_w_6k(matching_tests{ii,4}(kk),1);
                m = m_6k;
            end
            matching_tests{ii,7}(kk) = matching_tests{ii,6}(kk)/(f_w*d_sph); %Extracting reduced velocity
            kk = kk+1;
        end
    end
    
    for iii = [1,6]
        temp = find(matching_tests{ii,2}==iii);
        temp_nan = find(~(matching_tests{ii,2}==iii));
        uniq = unique(round(matching_tests{ii,7}(temp)/0.5)*0.5); %Finds all pump speeds that were tested
        uniq_idx = 1:length(uniq); %Number of different pump speeds
        for jjj = 1:length(uniq)
            samespring_ustar(temp) = matching_tests{ii,7}(temp);
            samespring_ustar(temp_nan) = NaN;
            temp2 = find(round(samespring_ustar/0.5)*0.5==uniq(jjj)); %Finds the indexes in the matching_tests that matches each unique pump frequency
            matching_tests{ii,5}(temp2) = uniq_idx(jjj); %Index in an array where each pump speed would belong
        end
    end
end

%% Data Processing
dt = 1/f_s;

[test_size ~] = size(matching_tests); %Gives the number of unique configurations that were tested
for ii=1:test_size
    for jj=1:length(matching_tests{ii}) %Gives the number of tests for each configuration
    f_pump = matching_tests{ii,3}(jj);
    U = matching_tests{ii,6}(jj);

    k_temp = matching_tests{ii,2}(jj);
    if k_temp == 1
        k = k_1k(matching_tests{ii,4}(jj));
        kk = 1;
        m_a = m_a_1k(matching_tests{ii,4}(jj));
        c = c_1k(matching_tests{ii,4}(jj));
        f_w = f_w_1k(matching_tests{ii,4}(jj),1);
        f_w_95 = f_w_1k_95;
    else
        k = k_6k(matching_tests{ii,4}(jj),1);
        kk = 2;
        m_a = m_a_6k(matching_tests{ii,4}(jj));
        c = c_6k(matching_tests{ii,4}(jj));
        f_w = f_w_6k(matching_tests{ii,4}(jj),1);
        f_w_95 = f_w_6k_95;
    end
    data = table2array(readtable(topfolder+matching_tests{ii,1}(jj)));
    time = data(:,1);
    encoder = data(:,2);
    encoder_offset = 0;
    if f_pump==0
        encoder_offset = mean(encoder);
    else
        encoder = encoder-encoder_offset;
    end

    iii = matching_tests{ii,4}(jj);
    jjj = matching_tests{ii,5}(jj);
    
    f_testing{ii,kk}(iii,jjj) = f_pump;
    U_testing{ii,kk}(iii,jjj) = U;
    u_red{ii,kk}(iii,jjj) = U/(f_w(1)*d_sph); %Indexing is {configuration, spring constant}(test number, pump speed)
    u_red_95{ii,kk}(iii,jjj) = sqrt((U_68/(f_w*d_sph))^2+(U*f_w_95/(2*f_w*d_sph))^2)*2;
    pump_f{ii,kk}(iii,jjj) = f_pump;

    clear data

    %% Filtering
    clear PSD PSD_filt PSD_SG
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
    clear f_peaks mx phase f_windowed A_y_max peak_idx PSD_freq
    % figure
    [mx,phase,f_windowed] = psdd3_sayre(f_s,encoder_filt,12000,6000,2);
    f = f_windowed;
    f_norm = f_windowed./f_w(1);
    meanpwr = mean(mx,2);

    if single_test == 1
        nfft = 500000;
        [PSD_freq, PSD_norm{ii,kk}(iii,jjj,:)] = norm_PSD_calc(f_s,encoder_filt,nfft,2);
        PSD_freq_norm{ii,kk}(iii,jjj,:) = PSD_freq/f_w(1);
    end

    [pwr_max peak_idx] = max(meanpwr); %Finds the max power and location after taking average)
    
    p_f = polyfit(f(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2);
    p_f_norm = polyfit(f_norm(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2); %Fits a polynomial to the top of the mean power
    f_peak_temp = -p_f(2)/(2*p_f(1));
    f_peak{ii,kk}(iii,jjj) = f_peak_temp; %Setting 1st derivative slope to be 0, finding the location of 0
    f_star_peak{ii,kk}(iii,jjj) = -p_f_norm(2)/(2*p_f_norm(1));

    if jjj==1
        f_star_peak{ii,kk}(iii,jjj) = NaN;
    end
    
    if ii==1 & iii==1
        f_vo_norm(iii) = (St*U/d_sph)/f_w(1);
    end
    
    u_norm{ii,kk}(iii,jjj) = (u_red{ii,kk}(iii,jjj)./f_star_peak{ii,kk}(iii,jjj))*St;
    u_norm_95{ii,kk}(iii,jjj) = sqrt((U_68*St/(f_peak_temp*d_sph))^2+(U*St_68/(f_peak_temp*d_sph))^2)*2;
    
    %% Data Processing
    velo = FivePointDiff(encoder_filt,f_s)';
    velo = filtfilt(b,a,velo);
    acc = FivePointDiff(velo,f_s)';
    acc = filtfilt(b,a,acc);
    A_rms = (rms(encoder_filt));
    y_max = max(abs(encoder_filt));
    
    pos_peaks = findpeaks(encoder_filt,'MinPeakProminence',0.01);
    neg_peaks = findpeaks(-encoder_filt,'MinPeakProminence',0.01);
    all_peaks = [pos_peaks; neg_peaks];
    peaks_percentile = prctile(all_peaks,[10 90]);
    
    data_length = length(acc);
    time = time(1:data_length);
    encoder_filt = encoder_filt(1:data_length);
    velo =  velo(1:data_length);

    F = m.*acc+c.*velo+k*encoder_filt; 
    F_pot = -m_a*acc;%-C_A*m_d*acc;
    F_vortex = F - F_pot;
    % F_vortex = (m+m_d).*acc+f_n(2).*velo+k*encoder_filt;
 
    force_norm = 0.5*rho*(U^2)*pi*d_sph^2/4;
    C_y = F/force_norm;
    C_pot = F_pot/force_norm;
    C_vortex = F_vortex/force_norm;
    
    % if kk == 1 && diagnose == true
    %     limits = 10000:100000;
    %     figure
    %     tiledlayout(4,1)
    %     nexttile
    %     plot(time(limits),encoder_filt(limits))
    %     ylabel('Displacement')
    %     title(u_red(ii,jj,kk,iii))
    %     nexttile
    %     plot(time(limits),C_y(limits))
    %     ylabel('C_y')
    %     nexttile
    %     plot(time(limits),C_vortex(limits))
    %     ylabel('C_{vortex}')
    %     nexttile
    %     plot(time(limits),C_pot(limits))
    %     ylabel('C_{pot}')
    % end

    C_y_rms{ii,kk}(iii,jjj) = rms(C_y-mean(C_y));
    C_pot_rms{ii,kk}(iii,jjj) = rms(C_pot-mean(C_pot));
    C_vortex_rms{ii,kk}(iii,jjj) = rms(C_vortex-mean(C_vortex));
    
    A_y_star{ii,kk}(iii,jjj) = sqrt(2)*A_rms/d_sph;
    peaks_10{ii,kk}(iii,jjj) = peaks_percentile(1)/d_sph;
    peaks_90{ii,kk}(iii,jjj) = peaks_percentile(2)/d_sph;
    pdicy{ii,kk}(iii,jjj) = sqrt(2)*A_rms./y_max; %Periodicity
    if jjj==1
        A_y_star{ii,kk}(iii,jjj) = NaN;
        pdicy{ii,kk}(iii,jjj) = NaN;
    end

    %% Phase Lag Calculations
    % encoder_hilbert = hilbert(encoder_filt);
    % C_y_hilbert = hilbert(C_y);
    % C_vortex_hilbert = hilbert(C_vortex);
    % % relative_phase_deg = 360/(2*pi)*(angle(C_y_hilbert)-angle(encoder_hilbert));
    % % relative_phase_deg = mod(relative_phase_deg, 360) - 180;
    % % relative_phase_deg = mod(relative_phase_deg,360);
    % % relative_phase_deg(relative_phase_deg>180) = 360- relative_phase_deg(relative_phase_deg > 180);
    % % figure
    % % tiledlayout(2,1)
    % % nexttile
    % % plot(relative_phase_deg)
    % C_y_phase(ii,jj,kk,iii) = mean(relative_phase_deg);
    % 
    % relative_phase_deg = 360/(2*pi)*(angle(C_vortex_hilbert)-angle(encoder_hilbert));
    % relative_phase_deg = mod(relative_phase_deg, 360) - 180;
    % % nexttile
    % % plot(relative_phase_deg)
    % C_vortex_phase(ii,jj,kk,iii) = mean(relative_phase_deg);

    [C_y_phase_alt{ii,kk}(iii,jjj), C_vortex_phase_alt{ii,kk}(iii,jjj)] = retrievephase2(f_s,f_peak{ii,kk}(iii,jjj),encoder_filt,C_y,C_vortex);
    [totalphase, vortexphase] = retrievephase1(encoder_filt,C_y,C_vortex);
    C_y_phase{ii,kk}(iii,jjj) = mean(totalphase);
    C_vortex_phase{ii,kk}(iii,jjj) = mean(vortexphase);
    if diagnose == true && kk==1
        figure
        plot(encoder_filt)
        hold on
        plot(C_vortex)
        ylabel('\phi')
        legend('Disp', 'Vortex Force')
        title(u_red{ii,kk}(iii,jjj))
    end
    end
% end
%% Determining the average and uncertainty bounds from the tests
if single_test == 1
    results = {u_red, u_norm, pdicy, C_y_rms, C_pot_rms, C_vortex_rms, C_y_phase, C_vortex_phase, A_y_star, f_star_peak, C_y_phase_alt, C_vortex_phase_alt, pump_f, u_red_95, u_norm_95, peaks_10, peaks_90, PSD_freq_norm, PSD_norm};
else
    results = {u_red, u_norm, pdicy, C_y_rms, C_pot_rms, C_vortex_rms, C_y_phase, C_vortex_phase, A_y_star, f_star_peak, C_y_phase_alt, C_vortex_phase_alt, pump_f, u_red_95, u_norm_95};
end

for kkk = 1:length(results)
    [results_ave{kkk}{ii}(:,:), results_upper{kkk}{ii}(:,:), results_lower{kkk}{ii}(:,:)] = ave_bounds(results,kkk,ii);
end

%% Plotting results from distances

%First plot is Ay_star

if single_test == 1
    % figure(freq_contour_fig)
    % plot_psd_fn(results_ave,1,16,17,ii,plot_legends,plotting_color)
    % if ii==1
    %     Ustar_temp = 0:23.5;
    %     f_vo_norm = St*Ustar_temp;
    %     plot(Ustar_temp,f_vo_norm,'k--','DisplayName','Static')
    %     yline(1,'k-','HandleVisibility','off')
    %     set(gca,'XMinorTick','on','YMinorTick','on')
    % end

    figure(A_y_star_pctile_fig)
    plot_fn_prc(results_ave,1,9,14,15,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style)
    if ii==1
        set(gca,'XMinorTick','on','YMinorTick','on')
        xlabel('$U^*$')
        ylabel('$A^*$')
        set(get(gca,'ylabel'),'rotation',0)
    end

end


figure(A_y_star_fig)
% plot_legends=0;
hold on
if ii==1
    if plot_reference == 1
        plot(u_red_A_star_sareen,A_star_sareen,'k-s','DisplayName','Sareen 2018b');
        plot(u_red_A_star_govwill,A_star_govwill,'k-d','DisplayName','Govhardan 2005');
    end
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('$U^*$')
    ylabel('$A^*$')
    set(get(gca,'ylabel'),'rotation',90)
end
plot_fn(results_ave,results_lower,results_upper,1,9,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,14)
% exportgraphics(A_y_star_fig,["figures\A_y_star_"+uniq_configs(ii)+".png"],'Resolution',300,'BackgroundColor', bgColor);
% set(get(gca,'ylabel'),'rotation',90)
% drawnow
% frame = getframe(gcf);
% im = frame2im(frame);
%     [imind, cm] = rgb2ind(im, 256);
% 
%     % Write to the GIF File
% if ii == 1
%     imwrite(imind, cm, ['figures\' 'A_y_star.gif'], 'gif', 'Loopcount', 0, 'DelayTime', 2);
% else
%     imwrite(imind, cm, ['figures\' 'A_y_star.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', 2);
% end

%Plots of normalized reduced velocity
figure(A_y_norm_fig)
hold on
if ii==1
    if plot_reference == 1
        plot(unorm_unorm_sareen,unorm_Astar_sareen,'k-s','DisplayName','Sareen 2018b');
        plot(unorm_unorm_govwill,unorm_Astar_govwill,'k-d','DisplayName','Govhardan 2005');
    end
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('$(U^*/f^*)St$')
    ylabel('$A^*$')
    set(get(gca,'ylabel'),'rotation',0)
end

plot_fn(results_ave,results_lower,results_upper,2,9,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,15)

% dim = [0.35 0.75 0.5 0.1];
% annotation('textbox',dim,'String','Mode II','FitBoxToText','on','EdgeColor','none','Interpreter','latex')
% dim = [0.17 0.5 0.5 0.1];
% annotation('textbox',dim,'String','Mode I','FitBoxToText','on','EdgeColor','none','Interpreter','latex')
% dim = [0.5 0.6 0.5 0.1];
% annotation('textbox',dim,'String','Mode III/Plateau','FitBoxToText','on','EdgeColor','none','Interpreter','latex')

%Plots of frequency ratio
figure(f_star_fig)
hold on
if ii==1
    Ustar_temp = 0:23.5;
    f_vo_norm = St*Ustar_temp;
    plot(Ustar_temp,f_vo_norm,'k-','DisplayName','Static')
    yline(1,'k--','HandleVisibility','off')
    set(gca,'XMinorTick','on','YMinorTick','on')
end

plot_fn(results_ave,results_lower,results_upper,1,10,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,14) %Freq Ratio plot
% set(gca)
xlabel('$U^*$')
ylabel('$f^*$')
ylim([0.9 1.2])
set(get(gca,'ylabel'),'rotation',0)
% exportgraphics(f_star_fig,["figures\f_star_"+uniq_configs(ii)+".png"],'Resolution',300,'BackgroundColor', bgColor)

%Plots of lift coefficient
hold on
if ii==1
    figure(total_force_fig)
    hold on
    if plot_reference==1
        plot(u_red_totalforce_govwill,totalforce_govwill,'k-d','DisplayName','Govhardan 2005');
    end
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('$U^*$')
    ylabel('$C^{\prime}_{total}$')
    set(get(gca,'ylabel'),'rotation',90)
    % legend

    figure(vortex_force_fig)
    hold on
    if plot_reference==1
        plot(u_red_vortexforce_govwill,vortexforce_govwill,'k-d','DisplayName','Govhardan 2005');
    end
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('$U^*$')
    ylabel('$C^{\prime}_{vortex}$')
    set(get(gca,'ylabel'),'rotation',90)
    % legend

    figure(total_phase_fig)
    hold on
    if plot_reference==1
        plot(u_red_totalphase_sareen,totalphase_sareen,'k-s','DisplayName','Sareen 2018b');
        plot(u_red_totalphase_govwill,totalphase_govwill,'k-d','DisplayName','Govhardan 2005');
    end
    set(gca,'XMinorTick','on')
    xlabel('$U^*$')
    ylabel('$\phi_{total}$')
    set(get(gca,'ylabel'),'rotation',90)
    yline(90,'k--','HandleVisibility','off')
    % legend

    figure(vortex_phase_fig)
    hold on
    if plot_reference == 1
        plot(u_red_vortexphase_sareen,vortexphase_sareen,'k-s','DisplayName','Sareen 2018b');
        plot(u_red_vortexphase_govwill,vortexphase_govwill,'k-d','DisplayName','Govhardan 2005');
    end
    set(gca,'XMinorTick','on')
    xlabel('$U^*$')
    ylabel('$\phi_{vortex}$')
    set(get(gca,'ylabel'),'rotation',90)
    yline(90,'k--','HandleVisibility','off')
    % legend
end

figure(total_force_fig)
hold on
plot_fn(results_ave,results_lower,results_upper,1,4,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,14)

figure(vortex_force_fig)
hold on
plot_fn(results_ave,results_lower,results_upper,1,6,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,14)

figure(total_phase_fig)
hold on
plot_fn(results_ave,results_lower,results_upper,1,7,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,14)
ylim([0 180])
yticks(0:15:180);  % Set ticks every 15 units
yticklabels({'', '', '', '', '', '', '90', '', '', '', '', '', '180'});  % Set labels only at 90 and 180
ax = gca;
% if ii==1
%     % errorbar(squeeze(results_ave{1}(ii,jj,:)),squeeze(results_ave{11}(ii,jj,:)),squeeze(results_lower{11}(ii,jj,:)),squeeze(results_upper{11}(ii,jj,:)),'k-s','MarkerFaceColor','g','DisplayName','Cross')
%     modeII_line(ii,jj) = axis_norm(squeeze(results_ave{1}(ii,jj,:)),squeeze(results_ave{7}(ii,jj,:)),90,xlimits_forces(1),xlimits_forces(2),ax);
% end

figure(vortex_phase_fig)
hold on
plot_fn(results_ave,results_lower,results_upper,1,8,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,14)
ylim([0 180])
yticks(0:15:180);  % Set ticks every 15 units
yticklabels({'', '', '', '', '', '', '90', '', '', '', '', '', '180'});  % Set labels only at 90 and 180
% errorbar(squeeze(results_ave{1}(ii,jj,:)),squeeze(results_ave{12}(ii,jj,:)),squeeze(results_lower{12}(ii,jj,:)),squeeze(results_upper{12}(ii,jj,:)),'k-s','MarkerFaceColor','g','DisplayName','Cross')

%Periodicity Plot
figure(pdicy_fig)
% plot_legends = 1;
plot_fn(results_ave,results_lower,results_upper,1,3,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,14)

xlabel('$U^*$')
ylabel('$P$')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'XMinorTick','on','YMinorTick','on')
ylim([0.4 1])
% clear results_upper results_lower results_ave results
end

figure(griffin_fig)
xscale("log");
scatter(griffin_massdamp_samp,griffin_Astar_samp,'kd','DisplayName','Govhardan 2005');
plot(griffin_massdamp_fit,griffin_Astar_fit,'k-','DisplayName','Govhardan 2005');
plot(mass_damp_1k,max(squeeze(results_ave{9}{1}(1,:,:))),'ko','MarkerFaceColor','k');
set(gca,'XMinorTick','on','YMinorTick','on')
xlabel('$(m^*+C_A)\zeta$')
ylabel('$A^*$')
set(get(gca,'ylabel'),'rotation',0)

%% Annotation Playing
% figure(A_y_star_fig)
% delete(findall(gcf,'type','annotation'))
% % xlim([2 7])
% ylim([0 1.0])
% dim = [0.35 0.75 0.5 0.1];
% annotation('textbox',dim,'String','Mode II','FitBoxToText','on','EdgeColor','none','Interpreter','latex')
% dim = [0.17 0.4 0.5 0.1];
% annotation('textbox',dim,'String','Mode I','FitBoxToText','on','EdgeColor','none','Interpreter','latex')
% dim = [0.6 0.6 0.5 0.1];
% annotation('textbox',dim,'String','Mode III/Plateau','FitBoxToText','on','EdgeColor','none','Interpreter','latex')
figure(total_force_fig)
legend('Location','southeast','NumColumns',3)
xlim([0 22.5])
ylim([-0.01 0.35])
%% Figure Saving
% figure(A_y_star_fig)
% delete(findall(gcf,'type','annotation'))
% arrow_x = [0.775 0.775];
% arrow_y = [0.85 0.35];
% arrow_anno = annotation('arrow', arrow_x, arrow_y);
% arrow_anno.Color = 'black';
% arrow_anno.LineWidth = 2;
% 
% figure(pdicy_fig)
% delete(findall(gcf,'type','annotation'))
% arrow_x = [0.775 0.775];
% arrow_y = [0.6 0.85];
% arrow_anno = annotation('arrow', arrow_x, arrow_y);
% arrow_anno.Color = 'black';
% arrow_anno.LineWidth = 2;

saveas(A_y_star_fig,['figures\' 'A_y_star.eps'],'epsc')
saveas(total_force_fig,['figures\' 'total_force.eps'],'epsc')
saveas(vortex_force_fig,['figures\' 'vortex_force.eps'],'epsc')
saveas(total_phase_fig,['figures\' 'total_phase.eps'],'epsc')
saveas(vortex_phase_fig,['figures\' 'vortex_phase.eps'],'epsc')
saveas(A_y_norm_fig,['figures\' 'A_y_norm.eps'],'epsc')
saveas(griffin_fig,['figures\' 'griffin.eps'],'epsc')
saveas(pdicy_fig,['figures\' 'pdicy.eps'],'epsc')
saveas(f_star_fig,['figures\' 'fstar.eps'],'epsc')
if single_test == 1
    saveas(A_y_star_pctile_fig,'A_y_star_pctile.eps','epsc')
    saveas(freq_contour_fig,'freq_contour.eps','epsc')
end

exportgraphics(A_y_star_fig,['figures\' 'A_y_star.png'],'Resolution',300,'BackgroundColor', bgColor);
exportgraphics(total_force_fig,['figures\' 'total_force.png'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(vortex_force_fig,['figures\' 'vortex_force.png'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(total_phase_fig,['figures\' 'total_phase.png'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(vortex_phase_fig,['figures\' 'vortex_phase.png'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(A_y_norm_fig,['figures\' 'A_y_norm.png'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(griffin_fig,['figures\' 'griffin.png'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(pdicy_fig,['figures\' 'pdicy.png'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(f_star_fig,['figures\' 'f_star.png'],'Resolution',300,'BackgroundColor', bgColor)
if single_test == 1
    exportgraphics(A_y_star_pctile_fig,['figures\' 'A_y_star_pctile.png'],'Resolution',300,'BackgroundColor', bgColor)
    exportgraphics(freq_contour_fig,['figures\' 'freq_contour.png'],'Resolution',300,'BackgroundColor', bgColor)
end
%% Testing force subfigure
% Copy the contents of each existing figure to the new combined figure
% Subplot 1
if single_test == 1
    phase_subplot_fig = figure;
    phase_subplot_fig.Position = [100 100 600 800];
    hold on;
    tl = tiledlayout(3,1);
    nexttile
    ax1 = get(A_y_star_fig, 'CurrentAxes');
    copyobj(allchild(ax1), gca);
    title(ax1.Title.String); % Copy title from original axes
    ylabel('$A^*$')
    xlabel('')
    xticklabels({})
    xlim(ax1.XLim)
    text(-0.15, 1.0, 'a)', 'Units', 'normalized', 'FontWeight', 'bold');
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(get(gca,'ylabel'),'rotation',0)
    box on
    % legend
    
    % Subplot 2
    nexttile
    ax2 = get(total_phase_fig, 'CurrentAxes');
    copyobj(allchild(ax2), gca);
    title(ax2.Title.String); % Copy title from original axes
    set(gca,'XMinorTick','on')
    ylabel('$\phi_{total}$')
    set(get(gca,'ylabel'),'rotation',90)
    ylim([0 200])
    xlim(ax1.XLim)
    xticklabels({})
    xlabel('')
    yticks(0:90:180);  % Set ticks every 15 units
    yticklabels({'','90', '180'});  % Set labels only at 90 and 180
    yline(90,'k--')
    text(-0.15, 1.0, 'b)', 'Units', 'normalized', 'FontWeight', 'bold');
    box on
    % legend
    
    % Subplot 3
    nexttile
    ax3 = get(vortex_phase_fig, 'CurrentAxes');
    copyobj(allchild(ax3), gca);
    title(ax3.Title.String); % Copy title from original axes
    set(gca,'XMinorTick','on')
    xlabel('$U/f_{n,w}D$')
    ylabel('$\phi_{vortex}$')
    set(get(gca,'ylabel'),'rotation',90)
    ylim([0 200])
    xlim(ax1.XLim)
    yticks(0:90:180);  % Set ticks every 15 units
    yticklabels({'', '90', '180'});  % Set labels only at 90 and 180
    yline(90,'k--')
    text(-0.15, 1.0, 'c)', 'Units', 'normalized', 'FontWeight', 'bold');
    drawnow
    % 
    % hold on
    % x=0.5;
    box on
    modeII_line_gov = axis_norm(u_red_totalphase_govwill,totalphase_govwill,90,ax,tl);
    modeII_line_sareen = axis_norm(u_red_totalphase_sareen,totalphase_sareen,90,ax,tl);
    position_temp = tl.InnerPosition;
    norm_bottom = position_temp(2);
    norm_top = position_temp(2)+position_temp(4);
    for ii = 1:length(results_ave{1,1})
        modeII_line = axis_norm((squeeze(results_ave{1,1}{ii})),(squeeze(results_ave{1,7}{ii})),90,ax,tl);
        annotation(phase_subplot_fig, 'line', [modeII_line modeII_line], [norm_bottom norm_top], 'Color', plotting_color(ii,:), 'LineWidth', 1.5,'LineStyle',':');
    end
    annotation(phase_subplot_fig, 'line', [modeII_line_gov modeII_line_gov], [norm_bottom norm_top], 'Color', 'k', 'LineWidth', 1.5,'LineStyle','--');
    annotation(phase_subplot_fig, 'line', [modeII_line_sareen modeII_line_sareen], [norm_bottom norm_top], 'Color', 'k', 'LineWidth', 1.5,'LineStyle','-.');
    % annotation(phase_subplot_fig, 'line', [modeII_line_0 modeII_line_0], [norm_bottom norm_top], 'Color', plotting_color(1,:), 'LineWidth', 1.5,'LineStyle',':');
    % annotation(phase_subplot_fig, 'line', [modeII_line_current modeII_line_current], [norm_bottom norm_top], 'Color', plotting_color(2,:), 'LineWidth', 1.5,'LineStyle','--');
    % annotation(phase_subplot_fig, 'line', [modeII_line(1,1) modeII_line(1,1)], [0 1], 'Color', 'k', 'LineWidth', 1.5,'LineStyle','-');
    % legend
    saveas(phase_subplot_fig,['figures\' 'phase_amp_fig.eps'],'epsc')
    exportgraphics(phase_subplot_fig,['figures\' 'phase_amp_fig.png'],'Resolution',300,'BackgroundColor', bgColor)
end
% 

%% Extra Plots
% figure(A_y_star_fig)
% xlim([4.4,10.6])
% legend('Visible','off')
% errorbars = findall(gca, 'Type', 'ErrorBar'); % Get all lines in current axes
% lines = findall(gca,'Type','Line');
% 
% for i = 1:length(lines)
%     set(lines(i),'Visible','off')
% end
% 
% for i = 1:length(errorbars)
%     set(errorbars(i), 'Visible', 'on')
%     errorbars(i).LineStyle = 'none'; % Set line style to solid
%     % if i < 5
%     %     set(errorbars(i), 'Visible', 'off')
%     % end
% end