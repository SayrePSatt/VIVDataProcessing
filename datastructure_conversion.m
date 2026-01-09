close all
clear all
clc

test_distratios = ["000" "015" "020" "025" "030" "040" "050" "060" "070" "100"];% "020" "030"];
test_diaratios = ["_00" "_10"]; %"06" "08"];
test_spring = ["1k"];
load("pumpFit_freq2velo.mat");
%% Experiment Specification
% datafolder = "E:\vivscratch_complete\";
topfolder = "E:\vivscratch_complete\";
freedecayfolder = topfolder;

rho = 998;
d_sph = 0.0889;  %Diameter of Sphere
m_1k = 2.458347-(2/3)*0.0029;
m_6k = 2.458347;    %Oscillating Mass. 2.4295 for 90mm setup, 1.916 for 80mm setup
m_d = (4/3)*pi*(d_sph/2)^3*rho+rho*0.005^2*pi*d_sph/4; %Displaced mass
f_s = 1000;     %Sampling Frequency
C_A = 0.5;     %Added mass coefficient

temp_1k = table2array(readtable(freedecayfolder+"freeDecay/1k_08_18_2025/freedecay_1k_air.dat"));
f_n_1k = temp_1k(1,:);
f_n_1k_95 = temp_1k(2,1);
zeta_a_1k_95(1) = temp_1k(2,2);
temp_1k = table2array(readtable(freedecayfolder+"freeDecay/1k_08_18_2025/freedecay_1k_water.dat"));
f_w_1k_95 = temp_1k(2,1);
f_w_1k = temp_1k(1,:);
zeta_w_1k_95(1) = temp_1k(2,2);

temp_1k = table2array(readtable(freedecayfolder+"freeDecay/1k_08_18_2025/freedecay_1k_air.dat"));
f_n_1k(2,:) = temp_1k(1,:);
f_n_1k_95(2) = temp_1k(2,1);
zeta_a_1k_95(2) = temp_1k(2,2);
temp_1k = table2array(readtable(freedecayfolder+"freeDecay/1k_08_18_2025/freedecay_1k_water.dat"));
f_w_1k_95(2) = temp_1k(2,1);
f_w_1k(2,:) = temp_1k(1,:);
zeta_w_1k_95(2) = temp_1k(2,2);

f_n_1k(3,:)=f_n_1k(2,:);
f_n_1k(4,:)=f_n_1k(2,:);
f_w_1k(3,:)=f_w_1k(2,:);
f_w_1k(4,:)=f_w_1k(2,:);
f_n_1k_95(3) = f_n_1k_95(2);
f_n_1k_95(4) = f_n_1k_95(2);
f_w_1k_95(3) = f_w_1k_95(2);
f_w_1k_95(4) = f_w_1k_95(2);

zeta_a_1k = f_n_1k(:,2);
zeta_a_1k_95(3) = zeta_a_1k_95(2);
zeta_a_1k_95(4) = zeta_a_1k_95(2);

zeta_w_1k = f_w_1k(:,2);
zeta_w_1k_95(3) = zeta_w_1k_95(2);
zeta_w_1k_95(4) = zeta_w_1k_95(2);

% Base uncertainties
U_sigma68 = mdl.MSE;
m_sigma68 = 0.5e-3;

f_n_1k_sigma68 = f_n_1k_95/2;
f_w_1k_sigma68 = f_w_1k_95/2;

%% File reading

datafolder = topfolder+"aftertare\";
all_files = dir(datafolder);

for ii = 3:length(all_files)-1
    temp_config = all_files(ii).name;
    configs(ii-2) = convertCharsToStrings(temp_config(4:14));
    % distances =
end

uniq_configs = unique(configs);
uniq_configs = uniq_configs(contains(uniq_configs,test_diaratios) & contains(uniq_configs,test_distratios) & contains(uniq_configs,test_spring)); %Selects only the configurations selected for testing
uniq_configs = flip(uniq_configs);
uniq_configs = circshift(uniq_configs,1);

matching_tests = {};


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

%% Reading Data
[test_size ~] = size(matching_tests); %Gives the number of unique configurations that were tested
for ii=1:test_size
    for jj=1:length(matching_tests{ii}) %Gives the number of tests for each configuration
    f_pump = matching_tests{ii,3}(jj);
    U = matching_tests{ii,6}(jj);

    k_temp = matching_tests{ii,2}(jj);
    if k_temp == 1
        f_w = f_w_1k(matching_tests{ii,4}(jj),1);
        f_w_68 = f_w_1k_sigma68(matching_tests{ii,4}(jj));
        zeta_w = zeta_w_1k(matching_tests{ii,4}(jj));
        zeta_w_68 = zeta_w_1k_95(matching_tests{ii,4}(jj))/2;

        f_n = f_n_1k(matching_tests{ii,4}(jj),1);
        f_n_68 = f_n_1k_sigma68(matching_tests{ii,4}(jj));
        zeta_a = zeta_a_1k(matching_tests{ii,4}(jj));
        zeta_a_68 = zeta_a_1k_95(matching_tests{ii,4}(jj))/2;
    else
        f_w = f_w_6k(matching_tests{ii,4}(jj),1);
        f_w_68 = f_w_6k_sigma68(matching_tests{ii,4}(jj));
        zeta_w = zeta_w_6k(matching_tests{ii,4}(jj));
        zeta_w_68 = zeta_w_6k_95(matching_tests{ii,4}(jj))/2;

        f_n = f_n_6k(matching_tests{ii,4}(jj),1);
        f_n_68 = f_n_6k_sigma68(matching_tests{ii,4}(jj));
        zeta_a = zeta_a_6k(matching_tests{ii,4}(jj));
        zeta_a_68 = zeta_a_6k_95(matching_tests{ii,4}(jj))/2;
    end
    data_file = datafolder+matching_tests{ii,1}(jj);
    file_metadata = dir(data_file);
    data = table2array(readtable(data_file));
    time = data(:,1);
    encoder = data(:,2);
    encoder_offset = 0;
    if f_pump==0
        encoder_offset = mean(encoder);
    else
        encoder = encoder-encoder_offset;
    end

    red_velo = round(U/(f_w(1)*d_sph)*2)/2;
    red_velo_test = sprintf('%04.1f',red_velo);
    f_pump_test = sprintf('%05.2f',f_pump);
    %% Header Writing

    creation_date = datestr(file_metadata.date, 'mm-dd-yyyy_HH-MM');
    test_name = matching_tests{ii,1}(jj);
    test_name = extractBetween(test_name,4,24);

    new_filename = creation_date+"_"+test_name+red_velo_test+"_"+f_pump_test+".dat";

    header = {
    'Date',                datestr(file_metadata.date, 'mm/dd/yyyy');
    'Time',                datestr(file_metadata.date, 'HH:MM');
    'Test Name',           '1k VIV tests, old tests that have been reformatted to new data structure';
    'Lab Name',            'EFDL';
    'Test Collector',      'Sayre Satterwhite';
    'Air Temperature (C)', 'NA';
    'Water Temperature Temperature (C)', 'NA';
    'Isolated Runs?',      'No';
    'Direction of Tests:', 'Unknown?'
    };

    colHeader = {'Sphere Diameter (m)', 'Mass (kg)', 'f_n,w (Hz)', 'zeta,w', 'f_n,a (Hz)', 'zeta,a'};
    structure_data = [
    d_sph    m    f_w    zeta_w    f_n    zeta_a;
    0.0005   0.005 f_w_68 zeta_w_68 f_n_68 zeta_a_68
            ];

    fid = fopen(topfolder+"aftertare_newstructure\"+new_filename, 'w');

    % Write header
    for i = 1:size(header, 1)
        fprintf(fid, '%s\t%s\n', header{i,1}, header{i,2});
    end
    
    % Blank line after header (optional)
    fprintf(fid, '\n');
    
    % Write column headers
    fprintf(fid, '%s\t', colHeader{1:end-1});
    fprintf(fid, '%s\n', colHeader{end});
    
    % Write data rows
    for i = 1:size(structure_data,1)
        fprintf(fid, '%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\n', structure_data(i,:));
    end
    
    fprintf(fid, '\n');
    
    newColHeader = {'Time (s)', 'Encoder (m)', 'Analog (V)', 'TTL', 'Trigger Sync'};
    fprintf(fid, '%s\t', newColHeader{1:end-1});
    fprintf(fid, '%s\n', newColHeader{end});

    % **Second table data**
    for i = 1:size(time,1)
        fprintf(fid, '%0.7f\t%0.7f\t%0.7f\t%0.7f\t%0.7f\n', time(i), encoder(i), 0, 0, 0);
    end

    fclose(fid);

    clear data
    end
end