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

%% Plot Testing
datafolder = "D:\EFDL\vivscratch\";
topfolder = datafolder+"testData\";
all_files = dir(topfolder);


for ii = 3:length(all_files)
    temp_config = all_files(ii).name;
    configs(ii-2) = convertCharsToStrings(temp_config(4:11));
    % distances =
end

uniq_configs = unique(configs);

for ii = 1:length(uniq_configs)
    uniq_dist(ii) = extractBetween(uniq_configs(ii),1,3);
    uniq_dia(ii) = extractBetween(uniq_configs(ii),6,7);
end


