%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: ReducedFreq_Response.m
% Version: 1
% Date: 04/28/2025
% Author: Sayre Satterwhite (sayreps@umich.edu)
% Description: Takes displacement data for a sphere in VIV and determines
% relavant quantities, compares with historical results
% This is adapted from a previous version for use with new data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
%% Experiment Specification

rho = 998;
d_sph = 0.08;  %Diameter of Sphere
m = 1.3928;    %Oscillating Mass
m_d = (4/3)*pi*(d_sph/2)^3*rho; %Displaced mass
f_s = 100;     %Sampling Frequency
C_A = 0.5;     %Added mass coefficient
f_n = table2array(readtable('free_decay_90mm_4k/freedecay_air.dat'));
f_n = f_n(1,:);
f_w = table2array(readtable('free_decay_90mm_4k/freedecay_water.dat'));
f_w = f_w(1,:);
m_a = ((f_n(1)/f_w(1))^2-1)*m;
St = 0.19;
k = (f_n(1)*2*pi)^2*m;
c = f_n(2)*2*sqrt((m+m_d*C_A)*k);
mass_damp = (m/m_d+C_A)*f_n(2);

%% Data Import
files = dir('RawData_90mm\*.csv');

for i=1:length(files)
    f_pump(i) = str2double(files(i).name(11:15));
    if f_pump(i) == 0
        U(i) = 0.0;
    else
        U(i) = pf2u(f_pump(i));
    end
    filename = fullfile(files(i).folder,files(i).name);
    data = import_data(filename);
    time{i} = data(:,1);
    encoder{i} = data(:,2);
    if i==1
        encoder_offset = mean(encoder{i});
    else
        encoder{i} = encoder{i}-encoder_offset;
    end
end
dt = 1/1000;

first = 1;
last = length(files);

u_red = U/(f_w(1)*d_sph); %Reduced Velocity

clear data files
%% Investigating Noise, Butterworth, S-G Filter

clear PSD PSD_filt PSD_SG
close all
subplot(3,1,1)
for i=1:length(U)
    n = length(time{i});
    fhat = fft(encoder{i}, n); % Compute the fast fourier transform
    PSD = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)
    
    freq = 1/(dt*n)*(0:n); % Create x-qxis of frequencies in Hz
    L = 1:floor(n/2); %only plot the first half of freqs
    plot(freq(L), PSD(L), 'LineWidth', 3); hold on
    set(gca, 'FontSize',12)
end
xlim([0.4 0.55])
xlabel('Frequency (f)')
ylabel('PSD')
xline(f_w(1))
title('Unfiltered Encoder Data')

f_c = 2;   %Cutoff frequency
[b,a] = butter(4,f_c/(f_s/2),"low");

for i=1:length(U)
    encoder_filt{i} = filtfilt(b,a,encoder{i}/1000);
end

subplot(3,1,2)
for i=1:length(U)
    n = length(time{i});
    fhat = fft(encoder_filt{i}, n); % Compute the fast fourier transform
    PSD_filt{i} = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)
    
    freq_filt{i} = 1/(dt*n)*(0:n); % Create x-qxis of frequencies in Hz
    L = 1:floor(n/2); %only plot the first half of freqs
    plot(freq_filt{i}(L), PSD_filt{i}(L), 'LineWidth', 3); hold on
    set(gca, 'FontSize',12)
end
xlim([0.4 0.55])
xlabel('Frequency (f)')
ylabel('PSD')
xline(f_w(1))
title("4th Order Butterworth Filter on Encoder Data")

subplot(3,1,3) %Savitzky-Golay Filter currently overfit
for i=1:length(U)
    PSD_SG{i} = sgolayfilt(PSD_filt{i},5,29);
    plot(freq_filt{i}(L),PSD_SG{i}(L),'LineWidth',3); hold on
    set(gca, 'FontSize',12)
end
xlim([0.4 0.55])
xlabel('Frequency (f)')
ylabel('PSD')
xline(f_w(1))
title("Savitzgy-Golay Filter")

%% Freq Investigation
close all
clear f_peaks mx phase f_windowed A_y_max peak_idx
figure
for i=1:length(U)
    [mx{i},phase{i},f_windowed{i}] = psdd3_sayre(f_s,encoder_filt{i},2048,1024,2);
    f{i} = f_windowed{i};
    f_norm{i} = f_windowed{i}./f_w(1);
    meanpwr = mean(mx{i},2);

    [pwr_max peak_idx] = max(meanpwr); %Finds the max power and location after taking average)
    
    p_f = polyfit(f{i}(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2);
    p_f_norm = polyfit(f_norm{i}(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2); %Fits a polynomial to the top of the mean power

    f_peaks(i) = -p_f(2)/(2*p_f(1)); %Setting 1st derivative slope to be 0, finding the location of 0
    f_star_peaks(i) = -p_f_norm(2)/(2*p_f_norm(1));
    % 
    % f_peaks{i}(j) = f_windowed{i}(peak_idx{i}(j),j);
    % 
    % f_peaks{i} = mean(f_peaks{i});
    % A_y_max{i} = mean(A_y_max{i});

    plot(f_windowed{i},mx{i},'LineWidth',3); hold on
    set(gca, 'FontSize',12)
end
xlim([0.4 0.55])
xlabel('Frequency')
ylabel('A_{rms}')
xline(f_w(1))
title('Frequency Response using Windowing')

% f_peaks = cell2mat(f_peaks);
% A_y_max = cell2mat(A_y_max);
% 
% f_star = f_peaks/f_w(1);

f_vo = St*U/d_sph;

U_normalized = (u_red./f_star_peaks)*St;

%% Data Processing


for i=1:length(U)
    velo{i} = FivePointDiff(encoder_filt{i},f_s)';
    acc{i} = FivePointDiff(velo{i},f_s)';
    A_rms(i) = (rms(encoder_filt{i}));
    y_max(i) = max(abs(encoder_filt{i}));
    
    data_length = length(acc{i});
    encoder_filt{i} = encoder_filt{i}(1:data_length,:);
    velo{i} =  velo{i}(1:data_length,:);

    %Investigating Strangeness in Phase
    mass_force{i} = m*acc{i};
    damping_force{i} = c.*velo{i};
    spring_force{i} = k*encoder_filt{i};

    F{i} = m.*acc{i}+c.*velo{i}+k*encoder_filt{i};
    F_pot{i} = -C_A*m_d*acc{i};
    F_vortex{i} = F{i} - F_pot{i};
    % F_vortex{i} = (m+m_d).*acc{i}+f_n(2).*velo{i}+k*encoder_filt{i};

    force_norm(i) = 0.5*rho*(U(i)^2)*(pi*(d_sph^2/4));
    C_y{i} = F{i}/force_norm(i);
    C_pot{i} = F_pot{i}/force_norm(i);
    C_vortex{i} = F_vortex{i}/force_norm(i);

    C_y_rms(i) = rms(C_y{i});
    C_vortex_rms(i) = rms(C_vortex{i});
end

A_y_star = sqrt(2)*A_rms/d_sph;
P = sqrt(2)*A_rms./y_max;

%% Force Investigation
% A_y_star = sqrt(2)*A_rms/d_sph;
% 
% figure
% scatter(u_red,A_y_star)
% 
% tiledlayout(1,3)
% nexttile
% plot(F{2}(1000:1100))
% nexttile
% plot(F_vortex{2}(1000:1100))
% nexttile
% plot(C_y{2}(1000:1100))

%% Determining Phase Lag and Periodicity

for i=1:length(U)-1
    % [vortex_crosscorr, vortex_lags] = xcorr(encoder_filt{i},C_vortex{i},'coeff');
    % [total_crosscorr, total_lags] = xcorr(encoder_filt{i},C_y{i},'coeff');
    % 
    % [~, vortex_maxIdx] = max(vortex_crosscorr);
    % [~, total_maxIdx] = max(total_crosscorr);
    % 
    % vortex_lagSamples = vortex_lags(vortex_maxIdx);
    % total_lagSamples = total_lags(total_maxIdx);
    % 
    % vortex_lagSeconds = vortex_lagSamples/dt;
    % total_lagSeconds = total_lagSamples/dt;
    % 
    % vortex_phaseLag(i) = (vortex_lagSeconds * 360)/f_peaks(i);
    % total_phaseLag(i) = (total_lagSeconds * 360)/f_peaks(i);
    encoder_hilbert = hilbert(encoder_filt{i});
    C_y_hilbert = hilbert(C_y{i});
    C_vortex_hilbert = hilbert(C_vortex{i});
    relative_phase_deg = rad2deg(angle(C_y_hilbert)-angle(encoder_hilbert));
    relative_phase_deg = mod(relative_phase_deg + 180, 360) - 180;
    C_y_phase(i) = mean(relative_phase_deg);

    relative_phase_deg = rad2deg(angle(C_vortex_hilbert)-angle(encoder_hilbert));
    relative_phase_deg = mod(relative_phase_deg + 180, 360) - 180;
    C_vortex_phase(i) = mean(relative_phase_deg);
end
%% Plots
sareen = csvread('Sareen_FS_Ustar.csv');
u_red_sareen = sareen(:,1);
A_star_sareen = sareen(:,2);

govwill =csvread('GovWill_UStar.csv');
u_red_govwill = govwill(:,1);
A_star_govwill = govwill(:,2);

subplot(3,1,1)
plot(u_red(2:12),A_y_star(2:12),'k-s','MarkerFaceColor','k','DisplayName','Current Study')
hold on
plot(u_red_sareen,A_star_sareen,'b-p','MarkerFaceColor','b','DisplayName','Sareen 2020');
plot(u_red_govwill,A_star_govwill,'r-d','MarkerFaceColor','r','DisplayName','Govhardan 2005');
set(gca, 'FontSize',12)
xlabel('U^*')
ylabel('A_{rms}^*')

% text(13.75,.7,'k\downarrow Sareen 2020')
% text(9,.775,'k\leftarrow Current Study')
% text(9.38,.39,'k\leftarrow Govhardan 2005')

xline([u_red(2),u_red(12)],'k--')

xlim([0 20])

% Customize x-axis
xticks(0:1:20); % Ticks every 1 unit
x_labels = cell(1, length(0:1:20));
x_label_vals = 0:5:20;

for i = 1:length(x_label_vals)
    x_labels{1 + (x_label_vals(i))} = num2str(x_label_vals(i));
end

xticklabels(x_labels);

% Customize y-axis
yticks(0:0.05:1); % Ticks every 0.05 unit

% Initialize a cell array for y-axis labels
y_labels = cell(size(0:0.05:1));

% Set labels at 0.25, 0.5, and 0.75
y_label_vals = [0.25, 0.5, 0.75];
for i = 1:length(y_label_vals)
    index = y_label_vals(i) * 20 + 1;  % Convert to index
    y_labels{index} = num2str(y_label_vals(i));
end

yticklabels(y_labels);

% Add annotations with text arrows for each line
% ax = gca;
% xlimits = ax.XLim;
% ylimits = ax.YLim;
% 
% Calculate normalized positions for arrow heads
% Example positions normalized based on axes limits
% xCurrentStudyNorm = (u_red(10) - xlimits(1)) / (xlimits(2) - xlimits(1));
% yCurrentStudyNorm = (A_y_star(10) - ylimits(1) - 0.05) / (ylimits(2) - ylimits(1));
% annotation('textarrow', [xCurrentStudyNorm - 0.1, xCurrentStudyNorm], ...
%     [yCurrentStudyNorm, yCurrentStudyNorm], 'String', 'Current Study', 'Color', 'k', 'FontSize', 10);
% 
% xSareenNorm = (u_red_sareen(50) - xlimits(1)) / (xlimits(2) - xlimits(1));
% ySareenNorm = (A_star_sareen(50) - ylimits(1)) / (ylimits(2) - ylimits(1));
% annotation('textarrow', [xSareenNorm, xSareenNorm], ...
%     [ySareenNorm-0.1, ySareenNorm], 'String', 'Sareen 2020', 'Color', 'k', 'FontSize', 10);
% 
% xGovwillNorm = (u_red_govwill(5) - xlimits(1)) / (xlimits(2) - xlimits(1));
% yGovwillNorm = (A_star_govwill(5) - ylimits(1)) / (ylimits(2) - ylimits(1));
% annotation('textarrow', [xGovwillNorm + 0.25, xGovwillNorm + 0.1], ...
%     [yGovwillNorm + 0.05, yGovwillNorm + 0.15], 'String', 'Govhardan 2005', 'Color', 'k', 'FontSize', 10);


% Phase Plots
subplot(3,1,2)
plot(u_red(2:12),C_y_rms(2:12),'k-s','MarkerFaceColor','k','DisplayName','C_{total}')
hold on
plot(u_red(2:12),C_vortex_rms(2:12),'b-p','MarkerFaceColor','b','DisplayName','C_{vortex}');
set(gca, 'FontSize',12)
xlabel('U^*')
ylabel('C')
legend

% text(13.75,.7,'k\downarrow Sareen 2020')
% text(9,.775,'k\leftarrow Current Study')
% text(9.38,.39,'k\leftarrow Govhardan 2005')
xlim([6 11])

%Freq Response
subplot(3,1,3)
plot(u_red(2:12),f_star_peaks(2:12),'k-s','MarkerFaceColor','k','DisplayName','f^*')
hold on
plot(u_red(2:12),f_vo(2:12)./f_w(1),'b-p','MarkerFaceColor','b','DisplayName','f_{vo}')
xlim([6 11])
set(gca, 'FontSize',12)
xlabel('U^*')
ylabel('f/f_n')
legend

figure
plot(u_red(2:12),P(2:12),'k-s','MarkerFaceColor','k','DisplayName','f^*')
xlim([6 11])
set(gca, 'FontSize',12)
xlabel('U^*')
ylabel('Periodicity')

figure
subplot(2,1,1)
plot(u_red(2:12),C_y_phase(2:12),'b-s','MarkerFaceColor','k','DisplayName','f^*')
xlim([6 11])
set(gca, 'FontSize',12)
xlabel('U^*')
ylabel('Total Phase')

subplot(2,1,2)
plot(u_red(2:12),C_vortex_phase(2:12),'b-s','MarkerFaceColor','k','DisplayName','f^*')
xlim([6 11])
set(gca, 'FontSize',12)
xlabel('U^*')
ylabel('Vortex Phase')

%% Testing
% figure
% run = 12;
% plot(F{run}/max(F{run}))
% hold on
% plot(F_vortex{run}/max(F_vortex{run}))
% plot(F_pot{run}/max(F_pot{run}))
% plot(encoder_filt{run}/max(encoder_filt{run}(1000:10000)))
% legend('Total', 'Vortex', 'Potential', 'Disp')
% %Plot above is to rudementarily investigate the phase difference between
% %the different forces
% title("Forces and Displacement for U^*="+u_red(run))
% 
% figure
% test = hilbert(encoder_filt{run});
% test2 = hilbert(F_vortex{run});
% mag = abs(test);
% phase = atan2(imag(test),real(test));
% relative_phase_deg = rad2deg(angle(test)-angle(test2));
% relative_phase_deg = mod(relative_phase_deg + 180, 360) - 180;
% plot(relative_phase_deg)
% title("Instantaneous Phase Lag for U^*="+u_red(run))
% 
% figure
% plot(mass_force{run})
% hold on
% plot(damping_force{run})
% plot(spring_force{run})
% plot(F_pot{run})
% legend('Interial','Damping','Spring','Potential')
% title("Breakdown of Forces for U^*="+u_red(run))
% 
% figure
