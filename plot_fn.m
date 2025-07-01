function plot_fn(results,lower_bound,upper_bound,x_num,y_num,config_idx,name,plot_legend,plotting_color)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
legend_label = strcat(sprintf('%02.1f',name),'D');

errorbar(squeeze(results{x_num}{config_idx}(:)),squeeze(results{y_num}{config_idx}(:)), ...
    squeeze(lower_bound{y_num}{config_idx}(:)),squeeze(upper_bound{y_num}{config_idx}(:)), ...
    'k-o','MarkerFaceColor',plotting_color(config_idx,:),'DisplayName',legend_label)

xlim([0 max(squeeze(results{x_num}{config_idx}(:)))+1])

if plot_legend == 1
    legend
end

end