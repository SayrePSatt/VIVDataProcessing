function [ave,upper,lower] = ave_bounds(quantity)
%ave_bounds takes in the desired quantity and produces the upper/lower
%bound lengths for use with error bar plotting
%   Detailed explanation goes here

quantity = squeeze(quantity);
ave = mean(quantity,1);
upper = max(quantity,[],1)-ave;
lower = ave-min(quantity,[],1);
end