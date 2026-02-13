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
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

%% Options for plotting
plot_legends = 1; %0 to not plot legends, 1 to plot legends
plot_reference = 0; %0 to not plot references
plot_errors = 0; %0 to not plot errorbars
single_test = 0; %Use for plotting the spectrogram curves and mean peaks curve
squareaxis = 0;

all_distratios = ["000" "015" "020" "025" "030" "040" "050" "060" "070" "100"];

test_distratios = ["000" "015" "020" "040" "070" "100"];
test_diaratios = ["_00" "_10"]; %"06" "08"];
test_spring = ["1k" "6k"];

[~, colormask, ~] = intersect(all_distratios,test_distratios);

bgColor = [255 255 255]/255;
figure_size = [100 100 600 350];
tick_size = [0.03 0.012];
%% Experiment Specification
% datafolder = "E:\vivscratch_complete\";
topfolder = "E:\EFDL\viv_newstructure\aftertare_newstructure\";

rho = 998;
C_A = 0.5;     %Added mass coefficient
St = 0.19;
St_68 = 0.005;

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
A_y_star_fig.Position = figure_size;
set(gca,'YLim',[0 1.5]);
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
hold on;

A_y_norm_fig = figure;
A_y_norm_fig.Position = figure_size;
set(gca,'XLim',[0.6 1.6]);
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
hold on;

phase_subplot_fig = figure;
phase_subplot_fig.Position = figure_size;
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
hold on;

total_force_fig = figure;
total_force_fig.Position = figure_size;
set(gca,'YLim',[0 0.5]);
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
hold on;

vortex_force_fig = figure;
vortex_force_fig.Position = figure_size;
set(gca,'YLim',[0 0.5]);
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
hold on;

total_phase_fig = figure;
total_phase_fig.Position = figure_size;
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
hold on;

vortex_phase_fig = figure;
vortex_phase_fig.Position = figure_size;
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
hold on;

pdicy_fig = figure;
pdicy_fig.Position = figure_size;
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
hold on;

f_star_fig = figure;
f_star_fig.Position = figure_size;
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
hold on;

test_fig = figure;
test_fig.Position = figure_size;
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
hold on;

griffin_fig = figure;
griffin_fig.Position = figure_size;
set(gca,'TickLength',tick_size);
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
if squareaxis == 1
    axis square
end
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
    configs(ii-2) = convertCharsToStrings(temp_config(18:25)); %18:28 to also include spring configs
    % distances =
end

uniq_configs = unique(configs);
uniq_configs = uniq_configs(contains(uniq_configs,test_diaratios) & contains(uniq_configs,test_distratios));% & contains(uniq_configs,test_spring)); %Selects only the configurations selected for testing
% uniq_configs_spacing = unique(extractBefore(uniq_configs,9));
% for ii=uniq_configs_spacing
% 
% end
uniq_configs = flip(uniq_configs);
uniq_configs = circshift(uniq_configs,1);

matching_tests = {};

extended_colormap = lines(7);
additional_colors = [204 0 204; %Magenta
                     255 204 102; %Yellow
                     255 0 255]/255; %Magenta
extended_colormap = [extended_colormap; additional_colors];
plotting_color(1,:) = [0 0 0];
plotting_color(2:length(all_distratios),:) = extended_colormap(1:length(all_distratios)-1,:); %lines(length(uniq_configs)-1); %zeros(length(uniq_configs)-1);% This is plotting with different colors
plotting_color = plotting_color(colormask,:);
% plotting_color(length(uniq_configs),:) = [1 0 0];
% plotting_color(2:end,:) = flipud(plotting_color(2:end,:));

marker_style = ["o"; "square"; "diamond"; "^"; "v"; ">"; "<"; "pentagram"; "hexagram";"*"];
marker_style = marker_style(1:length(all_distratios));
marker_style = marker_style(colormask);
% marker_style = flipud(marker_style(1:length(uniq_configs)));
% marker_style = circshift(marker_style,1,1);

uniq_configs = sort(uniq_configs); %Reorganizes data to be in order of separation distance

for ii = 1:length(uniq_configs)
    red_velo_est = [];
    filematch = [];
    red_velo_est_temp= [];
    uniq_dist(ii) = extractBetween(uniq_configs(ii),1,3); %Extracting distance ratios
    uniq_dia(ii) = extractBetween(uniq_configs(ii),6,7); %Extracting diameter ratios
    kk = 1;
    for jj = 3:length(all_files)
        filename = all_files(jj).name;
        if contains(filename,uniq_configs(ii)) && contains(filename,test_spring) && endsWith(filename,'.dat')
            red_velo_est_temp = str2double(cell2mat(extractBetween(filename,39,42)));
            if red_velo_est_temp == 0
                continue
            else
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

%% Data Processing

testing = [];
[num_uniq_configs, ~, ~] = size(matching_tests); %Gives the number of unique configurations that were tested
for ii=1:num_uniq_configs %each configuration
    [num_spring_configs, ~] = size(matching_tests{ii});
    for jj=1:num_spring_configs %Spring Config for each configuration
        num_red_velo = length(matching_tests{ii,jj});
        for kk=1:num_red_velo
            clear pdicy f_star_peak u_red u_red_68 A_y_star C_y_rms C_y_rms_68 C_pot_rms C_vortex_rms C_vortex_rms_68 C_y_phase C_vortex_phase f_vo_norm u_red_norm zeropad peaks_10 peaks_90
            num_datapoints = length(matching_tests{ii,jj}{kk});
            for iii = 1:num_datapoints
                data_idx = matching_tests{ii,jj}{kk}(iii);
                filename = all_files(data_idx).name;
                testing = [testing string(filename)];
                metadata = table2array(readtable(topfolder+filename,'Range','A12:F13'));
                data = table2array(readtable(topfolder+filename,'NumHeaderLines',14)); %Imports one file with corresponding data
                if data(end,1) > 225
                    data = data(50000:end,:);
                end
                %% Extracting metadata and run specifications
                f_pump = str2num(cell2mat(extractBetween(filename,44,48)));
                [U U_68_temp] = predict(mdl,f_pump,Alpha=0.05);%/(1.117645);
                U_68 = U-U_68_temp(1);
                d_sph = metadata(:,1);
                m_d = (4/3)*pi*(d_sph(1)/2)^3*rho+rho*0.005^2*pi*d_sph(1)/4; %Displaced mass
                m = metadata(:,2);
                f_nw = metadata(:,3);
                f_na = metadata(:,5);
                zeta = metadata(:,6);
                u_red(iii) = U/(f_nw(1)*d_sph(1));

                m_d = (4/3)*pi*(d_sph(1)/2)^3*rho+rho*0.005^2*pi*d_sph(1)/4;
                m_star = m(1)/m_d;
                C_A = ((f_na(1)./f_nw(1)).^2-1)*m_star;
                m_a = C_A*m_d;
                omega_na = 2*pi*f_na(1);
                k = m*omega_na.^2; %5.375; %(f_n(1)*2*pi)^2*m
                c = 2*m*omega_na*zeta(1);
                mass_damp = (m_star+C_A)*zeta(1);
                            %Make sure to include the calculations for
                            %uncertainty
                %Error propagation
                k_68 = sqrt((4*pi^2*f_na(1).^2*m(2)).^2+(8*pi^2*f_na(1)*m(1).*f_na(2)).^2);
                c_68 = sqrt((2*sqrt(m(1)*k(1)).*zeta(2)).^2 ...
                   +(zeta(1).*sqrt(k(1)/m(1)).*m(2)).^2 ...
                   +(zeta(1).*sqrt(m(1)./k(1)).*k_68).^2 );
                C_A_68 = sqrt((2*f_na(1).*m_star.*f_na(2)./f_nw(1)).^2 ...
                     +(-2*f_na(1)*m_star.*f_nw(2)./f_nw(1).^3).^2 ...
                     +((f_na(1)./f_nw(1)-1).*m(2)).^2);

                u_red_68(iii) = sqrt((U_68/(f_nw(1)*d_sph(1)))^2+(U*f_nw(2)/(2*f_nw(1)*d_sph(1)))^2);
                m_a_68 = C_A_68*m_d;

                %% Data processing
                time = data(:,1);
                f_s = 1/(time(2)-time(1));
                dt = 1/f_s;
                encoder = data(:,2);

                clear data metadata

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
                
                %% Frequency investigation
                clear f_peaks mx phase f_windowed A_y_max peak_idx %PSD_freq PSD_norm PSD_freq_norm
                % figure
                [mx,phase,f_windowed] = psdd3_sayre(f_s,encoder_filt,12000,6000,2);
                f = f_windowed;
                f_norm = f_windowed./f_nw(1);
                meanpwr = mean(mx,2);
                if single_test == 1
                    nfft = 500000;
                    [PSD_freq, PSD_norm(:,iii)] = norm_PSD_calc(f_s,encoder_filt,nfft,3*f_nw(1));
                    PSD_freq_norm(:,iii) = PSD_freq/f_nw(1);
                end
            
                [pwr_max peak_idx] = max(meanpwr); %Finds the max power and location after taking average)
                
                p_f = polyfit(f(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2);
                p_f_norm = polyfit(f_norm(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2); %Fits a polynomial to the top of the mean power
                f_peak_temp = -p_f(2)/(2*p_f(1));
                f_peak = f_peak_temp; %Setting 1st derivative slope to be 0, finding the location of 0
                f_star_peak(iii) = -p_f_norm(2)/(2*p_f_norm(1));

                f_vo_norm(iii) = (St*U/d_sph(1))/f_nw(1);
                
                u_norm(iii) = (u_red(iii)./f_star_peak(iii))*St;
                u_norm_68(iii) = sqrt((U_68*St/(f_peak_temp*d_sph(1)))^2+(U*St_68/(f_peak_temp*d_sph(1)))^2);

                %% Data Processing
                velo = FivePointDiff(encoder_filt,f_s)';
                % rms(velo)
                % velo = filtfilt(b,a,velo);
                acc = FivePointDiff(velo,f_s)';
                % acc = filtfilt(b,a,acc); %perhaps reintroduce these
                % additional filters if the signal becomes noisy, but
                % currently and additional delay may unnessacarily be
                % introduced
                A_rms = (rms(encoder_filt));
                y_max = max(abs(encoder_filt));
                
                pos_peaks = findpeaks(encoder_filt,'MinPeakProminence',0.01);
                neg_peaks = findpeaks(-encoder_filt,'MinPeakProminence',0.01);
                all_peaks = [pos_peaks; neg_peaks];
                peaks_percentile = prctile(all_peaks,[20 80]);
                
                data_length = length(acc);
                time = time(1:data_length);
                encoder_filt = encoder_filt(1:data_length);
                velo =  velo(1:data_length);
            
                F = m(1)*acc+c(1).*velo+k(1)*encoder_filt;
                F_pot = -m_a*acc;%-C_A*m_d*acc;
                F_vortex = F - F_pot;
                F_68 = sqrt((acc.*m(2)).^2+(velo.*c_68).^2+(encoder_filt.*k_68).^2);
                F_vortex_68 = sqrt((acc.*m(2)).^2+(acc.*m_a_68).^2+(velo.*c_68).^2+(encoder_filt.*k_68).^2);
             
                force_norm = 0.5*rho*(U^2)*pi*d_sph(1)^2/4;
                C_y = F/force_norm;
                C_pot = F_pot/force_norm;
                C_vortex = F_vortex/force_norm;
            
                C_y_68 = sqrt((-F.*U_68./(rho*U^3*pi*d_sph(1)^2)).^2 ...
                                   +(F_68./(0.5*rho*U^2*pi*d_sph(1)^2)).^2);
                C_vortex_68 = sqrt((-F_vortex.*U_68./(rho*U^3*pi*d_sph(1)^2)).^2 ...
                                       +(F_vortex_68./(0.5*rho*U^2*pi*d_sph(1)^2)).^2);

                C_y_rms(iii) = rms(C_y-mean(C_y));
                C_pot_rms(iii) = rms(C_pot-mean(C_pot));
                C_vortex_rms(iii) = rms(C_vortex-mean(C_vortex));
            
                C_y_rms_68(iii)= rms(C_y_68);%sum(C_y.*C_y_68)./(2*rms(C_y-mean(C_y))*sqrt(data_length));
                C_vortex_rms_68(iii) = rms(C_vortex_68);%sum(C_vortex.*C_vortex_68)./(2*rms(C_vortex-mean(C_vortex))*sqrt(data_length));

                A_y_star(iii) = sqrt(2)*A_rms/d_sph(1);
                peaks_10(iii) = peaks_percentile(1)/d_sph(1);
                peaks_90(iii) = peaks_percentile(2)/d_sph(1);
                pdicy(iii) = sqrt(2)*A_rms./y_max; %Periodicity

                %% Phase Lag Calculations
            
                [C_y_phase_alt, C_vortex_phase_alt] = retrievephase2(f_s,f_peak,encoder_filt,C_y,C_vortex);
                [totalphase, vortexphase] = retrievephase1(encoder_filt,C_y,C_vortex);
                C_y_phase(iii) = mean(totalphase);
                C_vortex_phase(iii) = mean(vortexphase);
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
            clear results
            %% Determining the average and uncertainty bounds from the tests
            zeropad = zeros(size(pdicy));
            % zeropad_psd = zeros(size(PSD_freq_norm));
            results = {[u_red; u_red_68], [u_norm; u_norm_68], [pdicy; zeropad], [C_y_rms; C_y_rms_68], [C_pot_rms; zeropad], [C_vortex_rms; C_vortex_rms_68], [C_y_phase; zeropad], [C_vortex_phase; zeropad], [A_y_star; zeropad], [f_star_peak; zeropad], [peaks_10; zeropad], [peaks_90; zeropad]};

            if single_test==1
                % psd_results = {PSD_freq_norm, PSD_norm};
                PSD_freq_norm_ave(:,kk) = mean(PSD_freq_norm,2);
                PSD_norm_ave(:,kk) = mean(PSD_norm,2);
            end

            for kkk = 1:length(results)
                [results_ave{kkk}{ii,jj}(kk), results_upper{kkk}{ii,jj}(kk), results_lower{kkk}{ii,jj}(kk)]= ave_bounds_newstructure(results{kkk});
            end

        end
    end
    %% Plotting Results
    
    if single_test == 1
        figure(freq_contour_fig)
        plot_psd_fn_newstructure(results_ave,1,PSD_freq_norm_ave,PSD_norm_ave,ii,plot_legends,plotting_color)
        if ii==1
            Ustar_temp = 0:23.5;
            f_vo_norm = St*Ustar_temp;
            plot(Ustar_temp,f_vo_norm,'k--','DisplayName','Static')
            yline(1,'k-','HandleVisibility','off')
            set(gca,'XMinorTick','on','YMinorTick','on','Layer','top')
        end
    
        figure(A_y_star_pctile_fig)
        plot_fn_prc(results_ave,1,9,11,12,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style)
        if ii==1
            set(gca,'XMinorTick','on','YMinorTick','on')
            xlabel('$U^*$')
            ylabel('$A^*$')
            set(get(gca,'ylabel'),'rotation',0)
        end
    end  
    
    %Plotting reduced amplitude
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
        set(get(gca,'ylabel'),'rotation',0)
    end
    plot_fn(results_ave,results_lower,results_upper,1,9,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,1)


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
    
    plot_fn(results_ave,results_lower,results_upper,2,9,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,2)
    xticks([0.6:0.2:1.6])
    xlim([0.6 1.6])

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
    
    plot_fn(results_ave,results_lower,results_upper,1,10,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,1) %Freq Ratio plot
    % set(gca)
    xlabel('$U^*$')
    ylabel('$f^*$')
    ylim([0.9 1.2])
    set(get(gca,'ylabel'),'rotation',0)

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
    plot_fn(results_ave,results_lower,results_upper,1,4,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,1)
    
    figure(vortex_force_fig)
    hold on
    plot_fn(results_ave,results_lower,results_upper,1,6,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,1)
    
    figure(total_phase_fig)
    hold on
    plot_fn(results_ave,results_lower,results_upper,1,7,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,1)
    % plot_fn(results_ave,results_lower,results_upper,1,11,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,14)
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
    plot_fn(results_ave,results_lower,results_upper,1,8,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,1)
    % plot_fn(results_ave,results_lower,results_upper,1,12,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,14)
    ylim([0 180])
    yticks(0:15:180);  % Set ticks every 15 units
    yticklabels({'', '', '', '', '', '', '90', '', '', '', '', '', '180'});  % Set labels only at 90 and 180
    % errorbar(squeeze(results_ave{1}(ii,jj,:)),squeeze(results_ave{12}(ii,jj,:)),squeeze(results_lower{12}(ii,jj,:)),squeeze(results_upper{12}(ii,jj,:)),'k-s','MarkerFaceColor','g','DisplayName','Cross')

    %Periodicity Plot
    figure(pdicy_fig)
    % plot_legends = 1;
    plot_fn(results_ave,results_lower,results_upper,1,3,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,1)
    
    xlabel('$U^*$')
    ylabel('$P$')
    set(get(gca,'ylabel'),'rotation',0)
    set(gca,'XMinorTick','on','YMinorTick','on')
    ylim([0.4 1])
end

figure(griffin_fig)
xscale("log");
scatter(griffin_massdamp_samp,griffin_Astar_samp,'kd','DisplayName','Govhardan 2005');
plot(griffin_massdamp_fit,griffin_Astar_fit,'k-','DisplayName','Govhardan 2005');
plot(mass_damp,max(squeeze(results_ave{9}{1}(1,:,:))),'ko','MarkerFaceColor','k');
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
% figure(total_force_fig)
% legend('Location','southeast','NumColumns',3)
% xlim([0 22.5])
% ylim([-0.01 0.35])
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

% saveas(A_y_star_fig,['figures\' 'A_y_star.eps'],'epsc')
% saveas(total_force_fig,['figures\' 'total_force.eps'],'epsc')
% saveas(vortex_force_fig,['figures\' 'vortex_force.eps'],'epsc')
% saveas(total_phase_fig,['figures\' 'total_phase.eps'],'epsc')
% saveas(vortex_phase_fig,['figures\' 'vortex_phase.eps'],'epsc')
% saveas(A_y_norm_fig,['figures\' 'A_y_norm.eps'],'epsc')
% saveas(griffin_fig,['figures\' 'griffin.eps'],'epsc')
% saveas(pdicy_fig,['figures\' 'pdicy.eps'],'epsc')
% saveas(f_star_fig,['figures\' 'fstar.eps'],'epsc')
% 
% saveas(A_y_star_pctile_fig,['figures\' 'A_y_star_pctile.eps'],'epsc')
% saveas(freq_contour_fig,'freq_contour.eps','epsc')

exportgraphics(A_y_star_fig,['figures\' 'A_y_star.pdf'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(total_force_fig,['figures\' 'total_force.pdf'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(vortex_force_fig,['figures\' 'vortex_force.pdf'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(total_phase_fig,['figures\' 'total_phase.pdf'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(vortex_phase_fig,['figures\' 'vortex_phase.pdf'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(A_y_norm_fig,['figures\' 'A_y_norm.pdf'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(griffin_fig,['figures\' 'griffin.pdf'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(pdicy_fig,['figures\' 'pdicy.pdf'],'Resolution',300,'BackgroundColor', bgColor)
exportgraphics(f_star_fig,['figures\' 'f_star.pdf'],'Resolution',300,'BackgroundColor', bgColor)
if single_test == 1
    exportgraphics(A_y_star_pctile_fig,['figures\' 'A_y_star_pctile.pdf'],'Resolution',300,'BackgroundColor', bgColor)
    exportgraphics(freq_contour_fig,['figures\' 'freq_contour.pdf'],'Resolution',300,'BackgroundColor', bgColor)
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
    ylabel('$\sqrt{2}A_{rms}/D$')
    xlabel('')
    xticklabels({})
    xlim(ax1.XLim)
    text(-0.15, 1.0, 'a)', 'Units', 'normalized', 'FontWeight', 'bold');
    set(gca,'XMinorTick','on','YMinorTick','on')
    set(get(gca,'ylabel'),'rotation',90)
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
    xlabel('$U^*$')
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
    legend
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