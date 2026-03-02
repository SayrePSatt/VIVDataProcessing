clear all %#ok<CLALL>
close all
clc
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');

%% Options for plotting
plot_legends = 1; %0 to not plot legends, 1 to plot legends
plot_reference = 0; %0 to not plot references
plot_errors = 0; %0 to not plot errorbars
single_test = 1; %Use for plotting the spectrogram curves and mean peaks curve
squareaxis = 0;
freq_plots = 1;
interpolate_scale = 8;

all_distratios = ["000" "015" "040"]% "020" "025" "030" "040" "050" "060" "070" "100"];

test_distratios = all_distratios;%["000" "015" "040"];
test_diaratios = ["_00" "_10"]; %"06" "08"];
test_spring = ["1k" "6k"];
% if length(test_spring)>1
%     multi_sprn

freq_cutoff = 7;

linesz = 1;

figure_width_spectrogram = 600;
figure_height_spectrogram = 200;
red_velo_forplotting = ["18.5_"];

[~, colormask, ~] = intersect(all_distratios,test_distratios);

bgColor = [255 255 255]/255;
figure_size = [100 100 600 350];
tick_size = [0.03 0.012];
size_marker = 6;
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

%% Figures
freq_PSD_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold off;

time_series_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

freq_contour_fig = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

freq_contour_fig_C_y = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;

freq_contour_fig_C_v = figure;
set(gcf, 'color', bgColor);
set(gca, 'color', bgColor);
hold on;
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
uniq_configs = circshift(uniq_configs,1)

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
        if contains(filename,uniq_configs(ii)) && contains(filename,test_spring) && endsWith(filename,'.dat')% && contains(filename,red_velo_forplotting)
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
    clear results_ave PSD_freq_norm_ave PSD_norm_ave PSD_freq_norm_ave_C_y PSD_norm_ave_C_y PSD_freq_norm_ave_C_v PSD_norm_ave_C_v
    [num_spring_configs, ~] = size(matching_tests{ii});
    for jj=1:num_spring_configs %Spring Config for each configuration
        num_red_velo = length(matching_tests{ii,jj});
        for kk=1:num_red_velo
            clear pdicy f_star_peak u_red u_red_68 u_norm u_norm_68 A_y_star C_y_rms C_y_rms_68 C_pot_rms C_vortex_rms C_vortex_rms_68 C_y_phase C_vortex_phase f_vo_norm u_red_norm zeropad peaks_10 peaks_90 e_vortex
            num_datapoints = length(matching_tests{ii,jj}{kk});
            for iii = 2 %1:num_datapoints
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
                k_spring = str2num(cell2mat(extractBetween(filename,27,27)));
                [U U_68_temp] = predict(mdl,f_pump,Alpha=0.05);%/(1.117645);
                U_68 = U-U_68_temp(1);
                d_sph = metadata(:,1);
                m_d = (4/3)*pi*(d_sph(1)/2)^3*rho+rho*0.005^2*pi*d_sph(1)/4; %Displaced mass
                m = metadata(:,2);
                f_nw = metadata(:,3);
                f_na = metadata(:,5);
                zeta = metadata(:,6);
                u_red = U/(f_nw(1)*d_sph(1));

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

                u_red_68 = sqrt((U_68/(f_nw(1)*d_sph(1)))^2+(U*f_nw(2)/(2*f_nw(1)*d_sph(1)))^2);
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

                % edot_vortex = C_vortex.*velo/d_sph(1);
                % e_vortex = trapz(time*f_peak,edot_vortex)/(max(time*f_peak)-min(time*f_peak));
            
                C_y_68 = sqrt((-F.*U_68./(rho*U^3*pi*d_sph(1)^2)).^2 ...
                                   +(F_68./(0.5*rho*U^2*pi*d_sph(1)^2)).^2);
                C_vortex_68 = sqrt((-F_vortex.*U_68./(rho*U^3*pi*d_sph(1)^2)).^2 ...
                                       +(F_vortex_68./(0.5*rho*U^2*pi*d_sph(1)^2)).^2);

                C_y_rms = rms(C_y-mean(C_y));
                C_pot_rms = rms(C_pot-mean(C_pot));
                C_vortex_rms = rms(C_vortex-mean(C_vortex));
            
                C_y_rms_68= rms(C_y_68);%sum(C_y.*C_y_68)./(2*rms(C_y-mean(C_y))*sqrt(data_length));
                C_vortex_rms_68 = rms(C_vortex_68);%sum(C_vortex.*C_vortex_68)./(2*rms(C_vortex-mean(C_vortex))*sqrt(data_length));

                A_y_star = sqrt(2)*A_rms/d_sph(1);
                peaks_10 = peaks_percentile(1)/d_sph(1);
                peaks_90 = peaks_percentile(2)/d_sph(1);
                pdicy = sqrt(2)*A_rms./y_max; %Periodicity

                %% Phase Lag Calculations
            
                % [C_y_phase_alt, C_vortex_phase_alt] = retrievephase2(f_s,f_peak,encoder_filt,C_y,C_vortex);
                % [totalphase, vortexphase] = retrievephase1(encoder_filt,C_y,C_vortex);
                % C_y_phase(iii) = mean(totalphase);
                % C_vortex_phase(iii) = mean(vortexphase);
                % if diagnose == true && kk==1
                %     figure
                %     plot(encoder_filt)
                %     hold on
                %     plot(C_vortex)
                %     ylabel('\phi')
                %     legend('Disp', 'Vortex Force')
                %     title(u_red{ii,kk}(iii,jjj))
                % end

                %% Frequency investigation
                clear f_peaks mx phase f_windowed A_y_max peak_idx PSD_freq PSD_norm PSD_freq_norm PSD_y %PSD_freq PSD_norm PSD_freq_norm
                % figure
                

                [mx,phase,f_windowed] = psdd3_sayre(f_s,encoder_filt,12000,6000,2);
                f = f_windowed;
                f_norm = f_windowed./f_nw(1);
                meanpwr = mean(mx,2);                

                [pwr_max peak_idx] = max(meanpwr); %Finds the max power and location after taking average)

                p_f = polyfit(f(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2);
                p_f_norm = polyfit(f_norm(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2); %Fits a polynomial to the top of the mean power
                f_peak_temp = -p_f(2)/(2*p_f(1));
                f_peak = f_peak_temp; %Setting 1st derivative slope to be 0, finding the location of 0
                f_star_peak(iii) = -p_f_norm(2)/(2*p_f_norm(1));
                % f_peak2 = f_star_peak(iii)*f_nw(1);
                f_vo_norm(iii) = (St*U/d_sph(1))/f_nw(1);

                if freq_plots == 1
                    nfft = 5000000;
                    [PSD_freq, PSD_y(:,iii), PSD_norm(:,iii)] = norm_PSD_calc(f_s,encoder_filt/d_sph(1),nfft,freq_cutoff*f_nw(1));
                    PSD_freq_norm(:,iii) = PSD_freq/f_peak;
                end

                % u_norm(iii) = (u_red(iii)./f_star_peak(iii))*St;
                % u_norm_68(iii) = sqrt((U_68*St/(f_peak_temp*d_sph(1)))^2+(U*St_68/(f_peak_temp*d_sph(1)))^2);
                %% Frequency Investigation of forces
                clear f_peaks mx phase f_windowed A_y_max peak_idx PSD_freq_C_y PSD_norm_C_y PSD_freq_norm_C_y PSD_freq_C_v PSD_norm_C_v PSD_freq_norm_C_v PSD_C_y PSD_C_v %PSD_freq PSD_norm PSD_freq_norm
                % figure
                % [mx,phase,f_windowed] = psdd3_sayre(f_s,C_y,12000,6000,2);
                % f = f_windowed;
                % f_norm = f_windowed./f_nw(1);
                % meanpwr = mean(mx,2);
                if freq_plots == 1
                    nfft = 50000;
                    [PSD_freq_C_y, PSD_C_y(:,iii), PSD_norm_C_y(:,iii)] = norm_PSD_calc(f_s,C_y,nfft,freq_cutoff*f_nw(1));
                    PSD_freq_norm_C_y(:,iii) = PSD_freq_C_y/f_peak;

                    [PSD_freq_C_v, PSD_C_v(:,iii), PSD_norm_C_v(:,iii)] = norm_PSD_calc(f_s,C_vortex,nfft,freq_cutoff*f_nw(1));
                    PSD_freq_norm_C_v(:,iii) = PSD_freq_C_v/f_peak;
                end
                % 
                % [pwr_max peak_idx] = max(meanpwr); %Finds the max power and location after taking average)
                % 
                % p_f = polyfit(f(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2);
                % p_f_norm = polyfit(f_norm(peak_idx-1:peak_idx+1),meanpwr(peak_idx-1:peak_idx+1),2); %Fits a polynomial to the top of the mean power
                % f_peak_temp = -p_f(2)/(2*p_f(1));
                % f_peak = f_peak_temp; %Setting 1st derivative slope to be 0, finding the location of 0
                % f_star_peak(iii) = -p_f_norm(2)/(2*p_f_norm(1));
                % 
                % f_vo_norm(iii) = (St*U/d_sph(1))/f_nw(1);

                %% Frequency Plots
                distance = round(str2double(char(test_distratios(ii)))/10*2)/2;
                if distance == 0
                    L_star = 'Isolated, '
                else
                    L_star = ['$L^*=$' +num2str(distance) ', ']
                end
                
                u_red_round = round(u_red*2)/2;
                figure(freq_PSD_fig)
                box on
                cla
                freq_PSD_fig.Position = [100 100 500 300];
                semilogy(PSD_freq_norm(:,iii),PSD_y(:,iii),'LineWidth',linesz,'Color','k','LineStyle','-')
                hold on
                semilogy(PSD_freq_norm_C_y(:,iii),PSD_C_y(:,iii),'LineWidth',linesz,'Color','r','LineStyle','-')
                semilogy(PSD_freq_norm_C_v(:,iii),PSD_C_v(:,iii),'LineWidth',linesz,'Color','b','LineStyle','-.')
                xlim([0 6]);
                ylim([10^-7 10^2]);
                xline([1 3 5],'LineStyle','--','Color',[80 80 80]/256)
                hold off
                xlabel('$f^*$')
                ylabel('PSD')
                yticks(logspace(-6,2,5))
                title([L_star '$U^*=$' num2str(u_red_round)])
                % % semilogy(PSD_freq_norm(:,iii),PSD_norm(:,iii),'k')
                % hold on
                % semilogy(PSD_freq_norm_C_y(:,iii),PSD_norm_C_y(:,iii),'r')
                % semilogy(PSD_freq_norm_C_v(:,iii),PSD_norm_C_v(:,iii),'b')
                % xlabel('$f^*$')
                % ylabel('dB/Hz')
                % xlim([0 7])
                % title("$U^*=$"+u_red)
                exportgraphics(freq_PSD_fig,['figures\spectralAnalysis\' filename(18:42) '_PSD.pdf'],'Resolution',300,'BackgroundColor', bgColor);

                %% Time Series Plots
                figure(time_series_fig)
                cla
                box on
                time_series_fig.Position = [100 100 500 300];
                cycles = 5;
                offset = 10; %How many seconds to wait
                cycle_data = int16(f_s*offset:f_s*(offset+round(cycles/f_nw(1))));
                tau = linspace(0,cycles,length(cycle_data));                
                title([L_star '$U^*=$' num2str(u_red_round)])
                plot(tau,encoder_filt(cycle_data)/d_sph(1),'LineWidth',linesz,'Color','k','LineStyle','-')
                hold on
                plot(tau,C_y(cycle_data),'LineWidth',linesz,'Color','r','LineStyle','-')
                plot(tau,C_vortex(cycle_data),'LineWidth',linesz,'Color','b','LineStyle','-.')
                xlabel('$\tau$')
                ylabel('$y^*,C_y,C_v$')
                ylim([-1.5 1.5]);
                
               exportgraphics(time_series_fig,['figures\spectralAnalysis\' filename(18:42) '_timeseries.pdf'],'Resolution',300,'BackgroundColor', bgColor);
            end
            
            results_ave{1}{ii,jj}(kk) = {[u_red]};
            results_ave{2}{ii,jj}(kk) = {[k_spring]};
            results_ave{3}{ii,jj}(kk) = {[f_peak]};
            % for kkk = 1:length(results)
            %     results_ave{kkk}{ii,jj}(kk) = ave_bounds_newstructure(results{kkk});
            % end

            if single_test==1
                % psd_results = {PSD_freq_norm, PSD_norm};
                PSD_freq_norm_ave{kk} = PSD_freq_norm(:,iii);
                PSD_norm_ave{kk} = PSD_norm(:,iii);

                PSD_freq_norm_ave_C_y{kk} = PSD_freq_norm_C_y(:,iii);
                PSD_norm_ave_C_y{kk} = PSD_norm_C_y(:,iii);

                PSD_freq_norm_ave_C_v{kk} = PSD_freq_norm_C_v(:,iii);
                PSD_norm_ave_C_v{kk} = PSD_norm_C_v(:,iii);
            end
        end
    end

    if freq_plots == 1
        figure(freq_contour_fig)
        freq_contour_fig.Position = [100 100 figure_width_spectrogram figure_height_spectrogram];
        box on
        plot_psd_fn_newstructure(interpolate_scale,results_ave,1,PSD_freq_norm_ave,PSD_norm_ave,ii,plot_legends,plotting_color)
        hold on
        if ii==1
        Ustar_temp = 0:22.5;
        f_vo_plot = St*Ustar_temp;
        end
        plot(Ustar_temp,f_vo_plot,'k-','DisplayName','Static')
        yline(1,'k--','HandleVisibility','off')
        set(gca,'XMinorTick','on','YMinorTick','on','Layer','top')
        xlim([0 22.5])
        title([L_star '$y^*$'])  
        ylabel('$f^*_{y^*}$')
        yticks(1:2:5)
        set(get(gca,'ylabel'),'rotation',0)
        exportgraphics(freq_contour_fig,['figures\spectralAnalysis\' convertStringsToChars(test_distratios(ii)) 'D_y_spectrogram.pdf'],'Resolution',500,'BackgroundColor', bgColor);
        exportgraphics(freq_contour_fig,['figures\spectralAnalysis\' convertStringsToChars(test_distratios(ii)) 'D_y_spectrogram.png'],'Resolution',500,'BackgroundColor', bgColor);
        savefig(freq_contour_fig,['figures\matlab_figs\' convertStringsToChars(test_distratios(ii)) 'D_y_spectrogram.fig']);
        hold off

        figure(freq_contour_fig_C_y)
        freq_contour_fig_C_y.Position = [100 100 figure_width_spectrogram figure_height_spectrogram];
        box on
        plot_psd_fn_newstructure(interpolate_scale,results_ave,1,PSD_freq_norm_ave_C_y,PSD_norm_ave_C_y,ii,plot_legends,plotting_color)
        hold on
        plot(Ustar_temp,f_vo_plot,'k-','DisplayName','Static')
        xlim([0 22.5])
        hold off
        yline(1,'k--','HandleVisibility','off')
        set(gca,'XMinorTick','on','YMinorTick','on','Layer','top')
        title([L_star '$C_{y}$'])
        ylabel('$f^*_{C_y}$')
        set(get(gca,'ylabel'),'rotation',0)
        yticks(1:2:5)
        exportgraphics(freq_contour_fig_C_y,['figures\spectralAnalysis\' convertStringsToChars(test_distratios(ii)) 'D_totalLift_spectrogram.pdf'],'Resolution',500,'BackgroundColor', bgColor);
        exportgraphics(freq_contour_fig_C_y,['figures\spectralAnalysis\' convertStringsToChars(test_distratios(ii)) 'D_totalLift_spectrogram.png'],'Resolution',500,'BackgroundColor', bgColor);
        savefig(freq_contour_fig_C_y,['figures\matlab_figs\' convertStringsToChars(test_distratios(ii)) 'D_totalLift_spectrogram.fig']);

        figure(freq_contour_fig_C_v)
        freq_contour_fig_C_v.Position = [100 100 figure_width_spectrogram figure_height_spectrogram];
        box on
        plot_psd_fn_newstructure(interpolate_scale,results_ave,1,PSD_freq_norm_ave_C_v,PSD_norm_ave_C_v,ii,plot_legends,plotting_color)
        hold on
        plot(Ustar_temp,f_vo_plot,'k-','DisplayName','Static')
        xlim([0 22.5])
        hold off
        yline(1,'k--','HandleVisibility','off')
        set(gca,'XMinorTick','on','YMinorTick','on','Layer','top')
        set(get(gca,'ylabel'),'rotation',0)
        title([L_star '$C_v$'])
        ylabel('$f^*_{C_v}$')
        yticks(1:2:5)
        exportgraphics(freq_contour_fig_C_v,['figures\spectralAnalysis\' convertStringsToChars(test_distratios(ii)) 'D_vortexLift_spectrogram.pdf'],'Resolution',500,'BackgroundColor', bgColor);
        exportgraphics(freq_contour_fig_C_v,['figures\spectralAnalysis\' convertStringsToChars(test_distratios(ii)) 'D_vortexLift_spectrogram.png'],'Resolution',500,'BackgroundColor', bgColor);
        savefig(freq_contour_fig_C_v,['figures\matlab_figs\' convertStringsToChars(test_distratios(ii)) 'D_vortexLift_spectrogram.fig']);

        % figure(A_y_star_pctile_fig)
        % plot_fn_prc(results_ave,1,9,11,12,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style)
        % if ii==1
        %     set(gca,'XMinorTick','on','YMinorTick','on')
        %     xlabel('$U^*$')
        %     ylabel('$A^*$')
        %     set(get(gca,'ylabel'),'rotation',0)
        % end
    end 
end