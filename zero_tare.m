%% Zero Extraction
% This script finds the mean for all tests and subtracts the appropriate
% amount

datafolder = "D:\EFDL\vivscratch_3\";
topfolder = datafolder+"6k_slow\";

all_files = dir(topfolder);


for ii = 3:length(all_files)
    temp_config = all_files(ii).name;
    configs(ii-2) = convertCharsToStrings(temp_config(4:14));
    % distances =
end

uniq_configs = unique(configs);

for ii = 1:length(uniq_configs)
    uniq_dist(ii) = extractBetween(uniq_configs(ii),1,3); %Extracting distance ratios
    uniq_dia(ii) = extractBetween(uniq_configs(ii),6,7); %Extracting diameter ratios
    kk = 1;
    for jj = 3:length(all_files)
        filename = all_files(jj).name;
        matching_tests{ii}(kk) = convertCharsToStrings(filename);
        kk=kk+1;
    end
end

[test_size ~] = size(matching_tests);

for ii=1:test_size
    for jj=1:length(matching_tests{ii}) %Gives the number of tests for each configuration
    f_pump = str2double(extractBetween(matching_tests{ii}(jj),25,29));
    clear data
    data = table2array(readtable(topfolder+matching_tests{ii,1}(jj)));
    time = data(:,1);
    encoder = data(:,2);
    if f_pump==0
        encoder_offset = mean(encoder);
    else
        encoder = encoder-encoder_offset;
    end

    data = [time encoder];
    if f_pump ~= 0.00
        writematrix(data,datafolder+"testDataZeroed\"+matching_tests{ii}(jj))
    end
    end
end