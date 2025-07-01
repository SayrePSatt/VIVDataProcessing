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

%% Reading Data
datafolder = "D:\EFDL\vivscratch\";
topfolder = datafolder+"testData\";
all_files = dir(topfolder);


for ii = 3:length(all_files)
    temp_config = all_files(ii).name;
    configs(ii-2) = convertCharsToStrings(temp_config(4:11));
    % distances =
end

uniq_configs = unique(configs);
matching_tests = {};

for ii = 1:length(uniq_configs)
    uniq_dist(ii) = extractBetween(uniq_configs(ii),1,3);
    uniq_dia(ii) = extractBetween(uniq_configs(ii),6,7);
    kk = 1;
    for jj = 3:length(all_files)
        filename = all_files(jj).name;
        if contains(filename,uniq_configs(ii)) && endsWith(filename,'.csv')
            matching_tests{ii,1}(kk) = convertCharsToStrings(filename);
            
            matching_tests{ii,2}(kk) = str2double(extractBetween(matching_tests{ii,1}(kk),13,13)); %Extracting the spring constant
            matching_tests{ii,3}(kk) = str2double(extractBetween(matching_tests{ii,1}(kk),25,29)); %Extracting Pump Speed
            matching_tests{ii,4}(kk) = str2double(extractBetween(matching_tests{ii,1}(kk),1,2)); %Extracting the test number
            kk = kk+1;
        end
    end
    
    for iii = [1, 6]
        temp = find(matching_tests{ii,2}==iii);
        uniq = unique(matching_tests{ii,3}(temp));
        uniq_idx = 1:length(uniq);
        for jjj = 1:length(uniq)
            temp2 = find(matching_tests{ii,3}==uniq(jjj));
            matching_tests{ii,5}(temp2) = uniq_idx(jjj);
        end
    end
end
