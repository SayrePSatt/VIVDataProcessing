function plot_fn(results,lower_bound,upper_bound,x_num,y_num,config_idx,name,plot_legend,plotting_color,plotting_marker)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
distance = str2double(extractBetween(name,1,3))/10;
legend_label = distance+"D";

errorbar(squeeze(results{x_num}{config_idx}(:)),squeeze(results{y_num}{config_idx}(:)), ...
    squeeze(lower_bound{y_num}{config_idx}(:)),squeeze(upper_bound{y_num}{config_idx}(:)), ...
    plotting_marker{config_idx},'MarkerFaceColor',plotting_color(config_idx,:), 'MarkerEdgeColor', ...
    'k','DisplayName',legend_label,LineStyle='none',Color='k',MarkerSize=8)

xlim([0 max(squeeze(results{x_num}{config_idx}(:)))+1])

if plot_legend == 1
    legend
end

end