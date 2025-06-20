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

%% Options for plotting
plot_legends = 1; %0 to not plot legends, 1 to plot legends

%% Experiment Specification
datafolder = "D:\EFDL\vivscratch\";

rho = 998;
d_sph = 0.0889;  %Diameter of Sphere
m = 2.42947;    %Oscillating Mass. 2.4295 for 90mm setup, 1.916 for 80mm setup
m_d = (4/3)*pi*(d_sph/2)^3*rho+0.005^2*pi*d_sph/4; %Displaced mass
f_s = 1000;     %Sampling Frequency
C_A = 0.5;     %Added mass coefficient
f_n = table2array(readtable(datafolder+"freeDecay/1k_06_19_2025/freedecay_1k_air.dat"));
f_n = f_n(1,:);
f_w = table2array(readtable(datafolder+"freeDecay/1k_06_19_2025/freedecay_1k_water.dat"));
f_w = f_w(1,:);
m_a = ((f_n(1)/f_w(1))^2-1)*m; %test
St = 0.19;
omegana = 2*pi*f_n(1);
k = m*omegana^2; %5.375; %(f_n(1)*2*pi)^2*m
c = f_n(2)*2*sqrt((m)*k);
m_star = m/m_d;
mass_damp = (m_star+C_A)*f_n(2);
scruton = 2*m*f_n(2)/(rho*d_sph^2); 

diagnose = false;

markers = ['s' 'd' '*'];
%% Importing external data for comparison
sareen = csvread('sareen2018b_ampphase.csv',3);
sareen(sareen==0) = NaN;
u_red_A_star_sareen = sareen(:,1);
A_star_sareen = sareen(:,2);
u_red_vortexphase_sareen = sareen(:,3);
vortexphase_sareen = sareen(:,4);
u_red_totalphase_sareen = sareen(:,5);
totalphase_sareen = sareen(:,6);

govwill = csvread('gov_2005_amplift.csv',3);
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

unorm_reference = csvread('unorm_reference.csv',3);
unorm_reference(unorm_reference==0) = NaN;
unorm_unorm_sareen = unorm_reference(:,1);
unorm_Astar_sareen = unorm_reference(:,2);
unorm_unorm_govwill = rmmissing(unorm_reference(:,3));
unorm_Astar_govwill = rmmissing(unorm_reference(:,4));

griffin_plot = csvread('griffin_govwill.csv',3);
griffin_plot(griffin_plot==0) = NaN;
griffin_massdamp_samp = griffin_plot(:,1);
griffin_Astar_samp = griffin_plot(:,2);
griffin_massdamp_fit = rmmissing(griffin_plot(:,3));
griffin_Astar_fit = rmmissing(griffin_plot(:,4));
%% Setting up plots for later

A_y_star_fig = figure;
hold on;

A_y_norm_fig = figure;
hold on;

phase_subplot_fig = figure;
hold on;

total_force_fig = figure;
hold on;

vortex_force_fig = figure;
hold on;

total_phase_fig = figure;
hold on;

vortex_phase_fig = figure;
hold on;

pdicy_fig = figure;
hold on;

f_star_fig = figure;
hold on;

test_fig = figure;
hold on;

griffin_fig = figure;
hold on;
%% Setting up folder directories
topfolder = datafolder+"testData\";
all_files = dir(topfolder+"*_diameter");
all_dir = all_files([all_files(:).isdir]);
num_diameters = numel(all_dir);

for ii=1:num_diameters
distance_files = dir(topfolder+all_dir(ii).name+"\*_distance");
distance_names(ii,:) = cellfun(@(x)str2num(x(1:3)),{distance_files.name})/10;
distance_dir = distance_files([distance_files(:).isdir]);
num_distances = numel(distance_dir);
plotting_color = lines(num_distances);

for jj=1:num_distances
test_files = dir(topfolder+all_dir(ii).name+"\"+distance_dir(jj).name+"\*_test");
test_dir = test_files([test_files(:).isdir]);
num_tests = numel(test_dir);

for kk=1:num_tests

filepath = topfolder+all_dir(ii).name+"\"+distance_dir(jj).name+"\"+test_dir(kk).name+"\";
files = dir(filepath+"*.csv");

for iii=1:length(files)
    dt = 1/f_s;
    f_pump = str2double(files(iii).name(18:22));
    if f_pump == 0
        U = 0.0;
    else
        U = pumpSpeedCalculator(f_pump);
    end
    filename = fullfile(files(iii).folder,files(iii).name);
    data = table2array(readtable(filename));
    time = data(:,1);
    encoder = data(:,2);
    if iii==1
        encoder_offset = mean(encoder);
    else
        encoder = encoder-encoder_offset;
    end

    u_red(ii,jj,kk,iii) = U/(f_w(1)*d_sph);
    pump_f(ii,jj,kk,iii) = f_pump;

    clear data

    %% Filtering
    clear PSD PSD_filt PSD_SG
    % subplot(3,1,1)
    n = length(time);
    fhat = fft(encoder, n); % Compute the fast fourier transform
    PSD = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)
    
    freq = 1/(dt*n)*(0:n); % Create x-qxis of frequencies in Hz
    L = 1:floor(n/2); %only plot the first half of freqs
    % plot(freq(L), PSD(L), 'LineWidth', 3); hold on
    % set(gca, 'FontSizf_ne',12)
    % 
    % xlim([0.4 0.55])
    % xlabel('Frequency (f)')
    % ylabel('PSD')
    % xline(f_w(1))
    % title('Unfiltered Encoder Data')
    
    f_c = 2;   %Cutoff frequency
    [b,a] = butter(4,f_c/(f_s/2),"low");
    
    encoder_filt = filtfilt(b,a,encoder);
    %% Freq Investigation
    clear f_peaks mx phase f_windowed A_y_max peak_idx
    % figure
    [mx,phase,f_windowed] = psdd3_sayre(f_s,encoder_filt,12000,6000,2);
    f = f_windowed;
    f_norm = f_windowed./f_w(1);
    meanpwr = mean(mx,2);

    [pwr_max peak_idx] = max(meanpwr); %Finds the max power and location after taking average)
    
    p_f = polyfit(f(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2);
    p_f_norm = polyfit(f_norm(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2); %Fits a polynomial to the top of the mean power

    f_peak(ii,jj,kk,iii) = -p_f(2)/(2*p_f(1)); %Setting 1st derivative slope to be 0, finding the location of 0
    f_star_peak(ii,jj,kk,iii) = -p_f_norm(2)/(2*p_f_norm(1));

    if iii==1
        f_star_peak(ii,jj,kk,iii) = NaN;
    end
    
    if ii==1 & jj==1 & kk==1
        f_vo_norm(iii) = (St*U/d_sph)/f_w(1);
    end
    
    u_norm(ii,jj,kk,iii) = (u_red(ii,jj,kk,iii)./f_star_peak(ii,jj,kk,iii))*St;
    
    %% Data Processing
    velo = FivePointDiff(encoder_filt,f_s)';
    velo = filtfilt(b,a,velo);
    acc = FivePointDiff(velo,f_s)';
    acc = filtfilt(b,a,acc);
    A_rms = (rms(encoder_filt));
    y_max = max(abs(encoder_filt));
    
    data_length = length(acc);
    time = time(1:data_length);
    encoder_filt = encoder_filt(1:data_length);
    velo =  velo(1:data_length);

    %Investigating Strangeness in Phase
    % mass_force = m*acc;
    % damping_force = c.*velo;
    % spring_force = k*encoder_filt;
    % 
    F = m.*acc+c.*velo/5+k*encoder_filt;
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

    C_y_rms(ii,jj,kk,iii) = rms(C_y-mean(C_y));
    C_pot_rms(ii,jj,kk,iii) = rms(C_pot-mean(C_pot));
    C_vortex_rms(ii,jj,kk,iii) = rms(C_vortex-mean(C_vortex));
    
    A_y_star(ii,jj,kk,iii) = sqrt(2)*A_rms/d_sph;
    pdicy(ii,jj,kk,iii) = sqrt(2)*A_rms./y_max; %Periodicity
    if iii==1
        A_y_star(ii,jj,kk,iii) = NaN;
        pdicy(ii,jj,kk,iii) = NaN;
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

    [C_y_phase_alt(ii,jj,kk,iii), C_vortex_phase_alt(ii,jj,kk,iii)] = retrievephase2(f_s,f_peak(ii,jj,kk,iii),encoder_filt,C_y,C_vortex);
    [totalphase, vortexphase] = retrievephase1(encoder_filt,C_y,C_vortex);
    C_y_phase(ii,jj,kk,iii) = mean(totalphase);
    C_vortex_phase(ii,jj,kk,iii) = mean(vortexphase);
    if diagnose == true && kk==1
        figure
        plot(encoder_filt)
        hold on
        plot(C_vortex)
        ylabel('\phi')
        legend('Disp', 'Vortex Force')
        title(u_red(ii,jj,kk,iii))
    end
end
end
%% Determining the average and uncertainty bounds from the tests

results = {u_red, u_norm, pdicy, C_y_rms, C_pot_rms, C_vortex_rms, C_y_phase, C_vortex_phase, A_y_star, f_star_peak, C_y_phase_alt, C_vortex_phase_alt, pump_f};
for jjj = 1:length(results)
    [results_ave{jjj}(ii,jj,:), results_upper{jjj}(ii,jj,:), results_lower{jjj}(ii,jj,:)] = ave_bounds(results{jjj}(ii,jj,:,:));
end

%% Plotting results from distances
figure(griffin_fig)
if ii==1 && jj==1
    xscale("log");
    scatter(griffin_massdamp_samp,griffin_Astar_samp,'kd','DisplayName','Govhardan 2005');
    plot(griffin_massdamp_fit,griffin_Astar_fit,'k-','DisplayName','Govhardan 2005');
    plot(mass_damp,max(squeeze(results_ave{9}(ii,jj,:))),'ko','MarkerFaceColor','k');
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('$(m^*+C_A)\zeta$')
    ylabel('$A^*$')
    set(get(gca,'ylabel'),'rotation',0)
end

%First plot is Ay_star
figure(A_y_star_fig)
hold on
if ii==1 && jj==1
    plot(u_red_A_star_sareen,A_star_sareen,'k-s','DisplayName','Sareen 2018b');
    plot(u_red_A_star_govwill,A_star_govwill,'k-d','DisplayName','Govhardan 2005');
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('$U^*$')
    ylabel('$A^*$')
    set(get(gca,'ylabel'),'rotation',0)
end

plot_fn(results_ave,results_lower,results_upper,1,9,ii,jj,distance_names(ii,jj),plot_legends,plotting_color)

%Plots of normalized reduced velocity
figure(A_y_norm_fig)
hold on
if ii==1 && jj==1
    plot(unorm_unorm_sareen,unorm_Astar_sareen,'k-s','DisplayName','Sareen 2018b');
    plot(unorm_unorm_govwill,unorm_Astar_govwill,'k-d','DisplayName','Govhardan 2005');
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('$(U^*/f^*)St$')
    ylabel('$A^*$')
    set(get(gca,'ylabel'),'rotation',0)
end

plot_fn(results_ave,results_lower,results_upper,2,9,ii,jj,distance_names(ii,jj),plot_legends,plotting_color)

dim = [0.55 0.8 0.5 0.1];
annotation('textbox',dim,'String','Mode II','FitBoxToText','on','EdgeColor','none','Interpreter','latex')
dim = [0.32 0.5 0.5 0.1];
annotation('textbox',dim,'String','Mode I','FitBoxToText','on','EdgeColor','none','Interpreter','latex')
dim = [0.75 0.65 0.5 0.1];
annotation('textbox',dim,'String','Mode III','FitBoxToText','on','EdgeColor','none','Interpreter','latex')

%Plots of frequency ratio
figure(f_star_fig)
hold on
if ii==1 && jj==1
    plot(squeeze(results_ave{1}(ii,jj,:)),f_vo_norm,'k-s','DisplayName',distance_dir(jj).name)
    set(gca,'XMinorTick','on','YMinorTick','on')
end

plot_fn(results_ave,results_lower,results_upper,1,10,ii,jj,distance_names(ii,jj),plot_legends,plotting_color)
% set(gca)
xlabel('$U^*$')
ylabel('$f^*$')
yline(1,'k--')
set(get(gca,'ylabel'),'rotation',0)

%Plots of lift coefficient
hold on
if ii==1 && jj==1
    figure(total_force_fig)
    hold on
    plot(u_red_totalforce_govwill,totalforce_govwill,'k-d','DisplayName','Govhardan 2005');
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('$U^*$')
    ylabel('$C_{total}$')
    set(get(gca,'ylabel'),'rotation',0)
    % legend

    figure(vortex_force_fig)
    hold on
    plot(u_red_vortexforce_govwill,vortexforce_govwill,'k-d','DisplayName','Govhardan 2005');
    set(gca,'XMinorTick','on','YMinorTick','on')
    xlabel('$U^*$')
    ylabel('$C_{vortex}$')
    set(get(gca,'ylabel'),'rotation',0)
    % legend

    figure(total_phase_fig)
    hold on
    plot(u_red_totalphase_sareen,totalphase_sareen,'k-s','DisplayName','Sareen 2018b');
    plot(u_red_totalphase_govwill,totalphase_govwill,'k-d','DisplayName','Govhardan 2005');
    set(gca,'XMinorTick','on')
    xlabel('$U^*$')
    ylabel('$\phi_{total}$')
    set(get(gca,'ylabel'),'rotation',0)
    yline(90,'k--','DisplayName','')
    % legend

    figure(vortex_phase_fig)
    hold on
    plot(u_red_vortexphase_sareen,vortexphase_sareen,'k-s','DisplayName','Sareen 2018b');
    plot(u_red_vortexphase_govwill,vortexphase_govwill,'k-d','DisplayName','Govhardan 2005');
    set(gca,'XMinorTick','on')
    xlabel('$U^*$')
    ylabel('$\phi_{vortex}$')
    set(get(gca,'ylabel'),'rotation',0)
    yline(90,'k--','DisplayName','')
    % legend
end

figure(total_force_fig)
hold on
plot_fn(results_ave,results_lower,results_upper,1,4,ii,jj,distance_names(ii,jj),plot_legends,plotting_color)

figure(vortex_force_fig)
hold on
plot_fn(results_ave,results_lower,results_upper,1,6,ii,jj,distance_names(ii,jj),plot_legends,plotting_color)

figure(total_phase_fig)
hold on
plot_fn(results_ave,results_lower,results_upper,1,7,ii,jj,distance_names(ii,jj),plot_legends,plotting_color)
ylim([0 180])
yticks(0:15:180);  % Set ticks every 15 units
yticklabels({'', '', '', '', '', '', '90', '', '', '', '', '', '180'});  % Set labels only at 90 and 180
ax = gca;
% errorbar(squeeze(results_ave{1}(ii,jj,:)),squeeze(results_ave{11}(ii,jj,:)),squeeze(results_lower{11}(ii,jj,:)),squeeze(results_upper{11}(ii,jj,:)),'k-s','MarkerFaceColor','g','DisplayName','Cross')
% modeII_line(ii,jj) = axis_norm(squeeze(results_ave{1}(ii,jj,:)),squeeze(results_ave{7}(ii,jj,:)),90,xlimits_forces(1),xlimits_forces(2),ax);

figure(vortex_phase_fig)
hold on
plot_fn(results_ave,results_lower,results_upper,1,8,ii,jj,distance_names(ii,jj),plot_legends,plotting_color)
ylim([0 180])
yticks(0:15:180);  % Set ticks every 15 units
yticklabels({'', '', '', '', '', '', '90', '', '', '', '', '', '180'});  % Set labels only at 90 and 180
% errorbar(squeeze(results_ave{1}(ii,jj,:)),squeeze(results_ave{12}(ii,jj,:)),squeeze(results_lower{12}(ii,jj,:)),squeeze(results_upper{12}(ii,jj,:)),'k-s','MarkerFaceColor','g','DisplayName','Cross')

%Periodicity Plot
figure(pdicy_fig)
plot_fn(results_ave,results_lower,results_upper,1,3,ii,jj,distance_names(ii,jj),plot_legends,plotting_color)

xlabel('$U^*$')
ylabel('$P$')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'XMinorTick','on','YMinorTick','on')
ylim([0 1])
clear results results_upper results_lower results_ave

end
end
%% Figure Saving
saveas(A_y_star_fig,'A_y_star.eps')
saveas(total_force_fig,'total_force.eps')
saveas(vortex_force_fig,'vortex_force.eps')
saveas(total_phase_fig,'total_phase.eps')
saveas(vortex_phase_fig,'vortex_phase.eps')
saveas(A_y_norm_fig,'A_y_norm.eps')
saveas(griffin_fig,'griffin.eps')

exportgraphics(A_y_star_fig,'A_y_star.jpg','Resolution',300);
exportgraphics(total_force_fig,'total_force.jpg','Resolution',300)
exportgraphics(vortex_force_fig,'vortex_force.jpg','Resolution',300)
exportgraphics(total_phase_fig,'total_phase.jpg','Resolution',300)
exportgraphics(vortex_phase_fig,'vortex_phase.jpg','Resolution',300)
exportgraphics(A_y_norm_fig,'A_y_norm.jpg','Resolution',300)
exportgraphics(griffin_fig,'griffin.jpg','Resolution',300)

%% Testing force subfigure
% Copy the contents of each existing figure to the new combined figure
% Subplot 1
phase_subplot_fig = figure;
hold on;
subplot(3, 1, 1, 'Parent', phase_subplot_fig);
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
legend

% Subplot 2
subplot(3, 1, 2, 'Parent', phase_subplot_fig);
ax2 = get(total_phase_fig, 'CurrentAxes');
copyobj(allchild(ax2), gca);
title(ax2.Title.String); % Copy title from original axes
set(gca,'XMinorTick','on')
ylabel('$\phi_{total}$')
set(get(gca,'ylabel'),'rotation',0)
ylim([0 180])
xlim(ax1.XLim)
xticklabels({})
xlabel('')
yticks(0:15:180);  % Set ticks every 15 units
yticklabels({'', '', '', '', '', '', '90', '', '', '', '', '', '180'});  % Set labels only at 90 and 180
yline(90,'k--')
text(-0.15, 1.0, 'b)', 'Units', 'normalized', 'FontWeight', 'bold');
legend

% Subplot 3
subplot(3, 1, 3, 'Parent', phase_subplot_fig);
ax3 = get(vortex_phase_fig, 'CurrentAxes');
copyobj(allchild(ax3), gca);
title(ax3.Title.String); % Copy title from original axes
set(gca,'XMinorTick','on')
xlabel('$U^*$')
ylabel('$\phi_{vortex}$')
set(get(gca,'ylabel'),'rotation',0)
ylim([0 180])
xlim(ax1.XLim)
yticks(0:15:180);  % Set ticks every 15 units
yticklabels({'', '', '', '', '', '', '90', '', '', '', '', '', '180'});  % Set labels only at 90 and 180
yline(90,'k--')
text(-0.15, 1.0, 'c)', 'Units', 'normalized', 'FontWeight', 'bold');
% 
% hold on
% x=0.5;
modeII_line_gov = axis_norm(u_red_totalphase_govwill,totalphase_govwill,90,ax1);
modeII_line_sareen = axis_norm(u_red_totalphase_sareen,totalphase_sareen,90,ax1);
annotation(phase_subplot_fig, 'line', [modeII_line_gov modeII_line_gov], [0 1], 'Color', 'k', 'LineWidth', 1.5,'LineStyle','--');
annotation(phase_subplot_fig, 'line', [modeII_line_sareen modeII_line_sareen], [0 1], 'Color', 'k', 'LineWidth', 1.5,'LineStyle','-.');
% annotation(phase_subplot_fig, 'line', [modeII_line(1,1) modeII_line(1,1)], [0 1], 'Color', 'k', 'LineWidth', 1.5,'LineStyle','-');
legend
saveas(phase_subplot_fig,'phase_amp_fig.eps')
saveas(phase_subplot_fig,'phase_amp_fig.jpg')


