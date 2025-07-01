function [ave,upper,lower] = ave_bounds(result,quantity,config)
%ave_bounds takes in the desired quantity and produces the upper/lower
%bound lengths for use with error bar plotting
%   Detailed explanation goes here

% quantity = squeeze(quantity);
quantity_1k = result{quantity}{config,1}(:,:);
ave_1k = mean(quantity_1k,1);
upper_1k = max(quantity_1k,[],1)-ave_1k;
lower_1k = ave_1k-min(quantity_1k,[],1);

quantity_6k = result{quantity}{config,2}(:,:);
ave_6k = mean(quantity_6k,1);
upper_6k = max(quantity_6k,[],1)-ave_6k;
lower_6k = ave_6k-min(quantity_6k,[],1);

ave = [ave_6k ave_1k];
upper = [upper_6k upper_1k];
lower = [lower_6k lower_1k];
end