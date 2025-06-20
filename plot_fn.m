function plot_fn(results,lower_bound,upper_bound,x_num,y_num,diameter_index,distance_index,name,plot_legend,plotting_color)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
legend_label = strcat(sprintf('%02.1f',name),'D');

errorbar(squeeze(results{x_num}(diameter_index,distance_index,:)),squeeze(results{y_num}(diameter_index,distance_index,:)), ...
    squeeze(lower_bound{y_num}(diameter_index,distance_index,:)),squeeze(upper_bound{y_num}(diameter_index,distance_index,:)), ...
    'k-o','MarkerFaceColor',plotting_color(distance_index,:),'DisplayName',legend_label)

xlim([0 max(squeeze(results{x_num}(diameter_index,distance_index,:)))+1])

if plot_legend == 1
    legend
end

end