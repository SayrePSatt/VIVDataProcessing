function plot_fn(results,lower_bound,upper_bound,x_num,y_num,config_idx,name,plot_legend,plotting_color,plotting_marker,errorbars,xerr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
distance = str2double(extractBetween(name,1,3))/10;
if distance == 0
    legend_label = "Iso.";
elseif distance == 10.5
    legend_label = "+0.5% f_n"
elseif distance == 9.5
    legend_label = "-0.5% f_n"
else
    legend_label = num2str(distance);
end

x_data = squeeze(results{x_num}{config_idx}(:));
y_data = squeeze(results{y_num}{config_idx}(:));
lwr_bound = squeeze(lower_bound{y_num}{config_idx}(:));
up_bound = squeeze(upper_bound{y_num}{config_idx}(:));
lwr_x_err = squeeze(lower_bound{xerr}{config_idx}(:));
up_x_err = squeeze(upper_bound{xerr}{config_idx}(:));
% xlen = length(x_data)
% ylen = length(y_data)
% ylower_len = length(lwr_bound)
% yupper_len = length(up_bound)
% xerr_len = length(x_err)

if errorbars == 1
    
    errorbar(x_data,y_data,lwr_bound,up_bound,lwr_x_err,up_x_err,...
        plotting_marker{config_idx}, ...
        'MarkerFaceColor',plotting_color(config_idx,:), ...[0 142 255]/255,
        'MarkerEdgeColor','k', ...
        'DisplayName',legend_label, ...
        LineStyle='none', ...
        LineWidth=1.2, ...
        Color=plotting_color(config_idx,:), ...
        MarkerSize=6)
else
    plot(x_data,y_data, ...
        plotting_marker{config_idx}, ...
        'MarkerFaceColor','none', ...%plotting_color(config_idx,:), ...
        'MarkerEdgeColor',plotting_color(config_idx,:), ...
        'DisplayName',legend_label, ...
        LineStyle='-', ...
        LineWidth=1.2, ...
        Color='w', ...
        MarkerSize=6)
end

xlim([0 max(squeeze(results{x_num}{config_idx}(:)))+1])

if plot_legend == 1
    plot_leg = legend('Location','eastoutside','NumColumns',1,'Interpreter', 'latex');
    set(plot_leg,'Color','none','Box','off')
    plot_leg.Title.String = '$L^*$';
end

box on

% Set major ticks every 5
set(gca, 'XTick', 0:5:24)

% Set minor ticks every 1
set(gca, 'XMinorTick', 'on')
% set(gca, 'XMinorTick', 0:1:24)

end