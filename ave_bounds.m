function [ave,upper,lower] = ave_bounds(result,quantity,config)
%ave_bounds takes in the desired quantity and produces the upper/lower
%bound lengths for use with error bar plotting
%   Detailed explanation goes here

% quantity = squeeze(quantity);
quantity_1k = result{quantity}{config,1}(:,:,:);
quantity_1k(quantity_1k==0) = NaN;
ave_1k = squeeze(mean(quantity_1k,1,'omitnan'));
size(ave_1k)
upper_1k = squeeze(max(quantity_1k,[],1,'omitnan'))-ave_1k;
lower_1k = ave_1k-squeeze(min(quantity_1k,[],1,'omitnan'));

quantity_6k = result{quantity}{config,2}(:,:,:);
quantity_6k(quantity_6k==0) = NaN;
ave_6k = squeeze(mean(quantity_6k,1,'omitnan'));
size(ave_6k)
upper_6k = squeeze(max(quantity_6k,[],1,'omitnan'))-ave_6k;
lower_6k = ave_6k-squeeze(min(quantity_6k,[],1,'omitnan'));

[rows_1k, columns_1k] = size(ave_1k);
[rows_6k, columns_6k] = size(ave_6k);

if columns_6k == columns_1k
    ave_6k = ave_6k';
    ave_1k = ave_1k';
    upper_6k = upper_6k';
    upper_1k = upper_1k';
    lower_6k = lower_6k';
    lower_1k = lower_1k';
end

ave = [ave_6k, ave_1k];
upper = [upper_6k, upper_1k];
lower = [lower_6k, lower_1k];
end