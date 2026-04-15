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
freq_plots = 0;

all_distratios = ["000" "015" "020" "025" "030" "040" "050" "060" "070" "100"]; %Do not change here, this controls the plot colors

test_distratios = ["000" "015" "020" "040" "070" "100"]; %Here, choose the 
test_diaratios = ["000" "010"]; %"06" "08"];
test_spring = ["1k" "6k"];

freq_cutoff = 6;

[~, colormask, ~] = intersect(all_distratios,test_distratios);

bgColor = [255 255 255]/255;
figure_size = [100 100 600 350];
tick_size = [0.03 0.012];
size_marker = 6;

lift_fig = figure;
%% Experiment Specification
% datafolder = "E:\vivscratch_complete\";
topfolder = "E:\test\";

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
    configs(ii-2) = convertCharsToStrings(temp_config(16:24)); %going to change after correction of time format
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
    span_dist_est = [];
    re_est = [];
    filematch = [];
    span_spac_est_temp= [];
    uniq_dist(ii) = extractBetween(uniq_configs(ii),1,3); %Extracting distance ratios
    uniq_dia(ii) = extractBetween(uniq_configs(ii),6,8); %Extracting diameter ratios
    kk = 1;
    for jj = 3:length(all_files)
        filename = all_files(jj).name;
        if contains(filename,uniq_configs(ii)) && ~contains(filename,'ZERO') && endsWith(filename,'.dat')
            span_dist_est_temp = str2double(cell2mat(extractBetween(filename,26,29)));
            span_dir = cell2mat(extractBetween(filename,32,35)); %get starbord/port
            if span_dir == "Port"
                span_dist_est_temp = -span_dist_est_temp;
            end
            span_dist_est = [span_dist_est span_dist_est_temp];

            re_est_temp = str2double(cell2mat(extractBetween(filename,39,43)));
            re_est = [re_est re_est_temp];
            filematch = [filematch jj];
        end
    end
    uniq_span = unique(span_dist_est);
    uniq_re = unique(re_est);
    for jjj = 1:length(uniq_re)
        for kk = 1:length(uniq_span)
            uniq_span(kk)
            temp_idx = find(uniq_span(kk)==span_dist_est & uniq_re(jjj)==re_est);
            matching_tests{ii,jjj}{kk} = filematch(temp_idx); %matching test indexing is: {configuration, reynolds, streamwise span, matching indicies for each est. red. velo}
        end
    end
end

%% Data processing
testing = [];
[num_uniq_configs, num_re_configs, ~] = size(matching_tests); %Gives the number of unique configurations that were tested
for ii=1:num_uniq_configs %each configuration
    for jj=1:num_re_configs %reynolds for each for each configuration
        num_span = length(matching_tests{ii,jj});
        for kk=1:num_span
            num_datapoints = length(matching_tests{ii,jj}{kk});
            for iii = 1:num_datapoints
                clear C_d C_l span_spac_ratio zeropad
                data_idx = matching_tests{ii,jj}{kk}(iii);
                filename = all_files(data_idx).name
                testing = [testing string(filename)];
                metadata = table2array(readtable(topfolder+filename,'Range','A10:G10')); %Imports the test specifications
                data = table2array(readtable(topfolder+filename,'NumHeaderLines',11)); %Imports one file with corresponding data
                if data(end,1) > 225
                    data = data(50000:end,:);
                end
                %% Extracting metadata and run specifications
                f_pump = 3.00; %str2num(cell2mat(extractBetween(filename,44,48))); NEED TO FIX FROM NEW FILENAME
                [U U_68_temp] = predict(mdl,f_pump,Alpha=0.05);
                d_sph = metadata(:,1);
                d_sphUS = metadata(:,2);
                rod_d = metadata(:,3);
                str_spac_ratio = metadata(:,4);
                span_spac_ratio(iii) = metadata(:,5)
                % span_dir = 
                renolds = metadata(:,7);

                %% Data Processing
                time = data(:,1);
                f_s = 1/(time(2)-time(1));
                dt = 1/f_s;
                force_transducer = data(:,2:7);
                clear data metadata
                
                force_transducer_ave = mean(force_transducer,1);
                force_norm = 0.5*rho*(U^2)*pi*d_sph(1)^2/4;

                C_d(iii) = kk;%force_transducer_ave(1)*force_norm;
                C_l(iii) = kk;%force_transducer_ave(2)*force_norm;
            end
            clear results
            zeropad = zeros(size(C_d));
            results = {[span_spac_ratio; zeropad], [C_d; zeropad],[C_l; zeropad]};
            for kkk = 1:length(results)
                [results_ave{kkk}{ii,jj}(kk), results_upper{kkk}{ii,jj}(kk), results_lower{kkk}{ii,jj}(kk)]= ave_bounds_newstructure(results{kkk});
            end
            figure(lift_fig)
            % plot_legends=0;
            hold on

        end
        plot_fn(results_ave,results_lower,results_upper,1,2,ii,uniq_configs(ii),plot_legends,plotting_color,marker_style,plot_errors,1)
    end
end