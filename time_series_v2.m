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

test_distratios = all_distratios;%["000" "015" "020" "025" "040" "100"];
test_diaratios = ["_00" "_10"]; %"06" "08"];
test_spring = ["1k"];
test_redvelo = ["_08.5" "_14.5" "_20.5"];
use_datapoint = 2;

timerseries_cycles = 50;
lissajous_cycles = 50;
offset = 30;

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

% extended_colormap = lines(7);
% additional_colors = [204 0 204; %Magenta
%                      255 204 102; %Yellow
%                      255 0 255]/255; %Magenta
% extended_colormap = [extended_colormap; additional_colors];
% plotting_color(1,:) = [0 0 0];
% plotting_color(2:length(all_distratios),:) = extended_colormap(1:length(all_distratios)-1,:); %lines(length(uniq_configs)-1); %zeros(length(uniq_configs)-1);% This is plotting with different colors
% plotting_color = plotting_color(colormask,:);
% 
% marker_style = ["o"; "square"; "diamond"; "^"; "v"; ">"; "<"; "pentagram"; "hexagram";"*"];
% marker_style = marker_style(1:length(all_distratios));
% marker_style = marker_style(colormask);
% % marker_style = flipud(marker_style(1:length(uniq_configs)));
% % marker_style = circshift(marker_style,1,1);

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
        if contains(filename,uniq_configs(ii)) && contains(filename,test_spring) && contains(filename,test_redvelo) && endsWith(filename,'.dat')
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
    uniq_dist = extractBetween(uniq_configs(ii),1,3);
    for jj=1:num_spring_configs %Spring Config for each configuration
        num_red_velo = length(matching_tests{ii,jj});
        for kk=1:num_red_velo
            clear pdicy f_star_peak u_red u_red_68 A_y_star C_y_rms C_y_rms_68 C_pot_rms C_vortex_rms C_vortex_rms_68 C_y_phase C_vortex_phase f_vo_norm u_red_norm zeropad peaks_10 peaks_90
            num_datapoints = length(matching_tests{ii,jj}{kk});
            for iii = use_datapoint
                data_idx = matching_tests{ii,jj}{kk}(iii);
                filename = all_files(data_idx).name;
                run = convertCharsToStrings(filename);
                testing = [testing string(filename)];
                metadata = table2array(readtable(topfolder+filename,'Range','A12:F13'));
                data = table2array(readtable(topfolder+filename,'NumHeaderLines',14)); %Imports one file with corresponding data
                if data(end,1) > 225
                    data = data(50000:end,:);
                end
                %% Extracting metadata and run specifications
                f_pump = str2num(cell2mat(extractBetween(filename,44,48)));
                red_velo_est = str2num(cell2mat(extractBetween(filename,39,42)));
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
                
                f_peak = -p_f(2)/(2*p_f(1));
                T = 1/f_peak;
                
                idx = (offset*T-time(1))<time & time<((offset+timerseries_cycles)*T(1));
                t = time(idx);
                t = t-t(1);
                encoder_filt_cont = encoder_filt(idx)/d_sph(1);

                idx_lissajous = (offset*T-time(1))<time & time<((offset+lissajous_cycles)*T(1));
                t_lissajous = time(idx_lissajous);
                t_lissajous = t_lissajous-t_lissajous(1);
                encoder_filt_lissajous = encoder_filt(idx_lissajous)/d_sph(1);

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
                
                C_y_lissajous = C_y(idx_lissajous);

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

                %% Plotting Time Series
                close all
                timeSeries = figure;
                % subplot(3,1,1)
                set(gcf, 'color', bgColor);
                set(gca, 'color', bgColor);
                timeSeries.Position = [100 100 500 250];
                plot(t/T,encoder_filt_cont,'k-')
                xlim([0 timerseries_cycles])
                limits = ceil(max(abs(encoder_filt_cont))/0.5)*0.5;
                % limits = 0.1;
                ylim([-limits limits])
                yticks([-limits 0 limits])
                % xlabel('$t/T$')
                ylabel('$y/D$')
                xlabel('$t/T$')
                distance = str2double(char(uniq_dist))/10;
                if distance == 0
                    L_star = 'Isolated ';
                else
                    L_star = ['$L^*=$' +num2str(distance) ' '];
                end
                title([L_star ' $U^*=$' num2str(red_velo_est)])
                set(gca,'xticklabel',0:10:timerseries_cycles)
                set(gcf, 'color', bgColor);
                set(gca, 'color', bgColor);

    
                %% Plotting Lissajous Curves
                lissajous_fig = figure;
                
                x = encoder_filt_lissajous;
                y = C_y_lissajous;
                
                plot(x,y,'LineWidth',1.5,'Color','k')
                xlabel('$y/D$')
                ylabel('$C_y$')
                title([L_star ' $U^*=$' num2str(red_velo_est)])
                xbounds = max(x+0.1);
                ybounds = max(y+0.1);
                xlim([-xbounds xbounds]);
                ylim([-ybounds ybounds]);
                axis square
                set(gcf, 'color', bgColor);
                set(gca, 'color', bgColor);
                figsize = get(gcf,'Position');
    
                %% Figure Export
                exportgraphics(timeSeries,['figures\timeSeries\' filename '_timeSeries.pdf'],'Resolution',300,'BackgroundColor', bgColor);
                exportgraphics(timeSeries,['figures\lissajous\' filename '_lissajous.pdf'],'Resolution',300,'BackgroundColor', bgColor);
            end
            clear results
        end
    end
end