function plot_fn_prc(results,x_num,y_num,prc_10,prc_90,config_idx,name,plot_legend,plotting_color,plotting_marker)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

distance = str2double(extractBetween(name,1,3))/10;
legend_label = distance+"D";

x_data = squeeze(results{x_num}{config_idx}(:));
y_data = squeeze(results{y_num}{config_idx}(:));
prc_10_data = squeeze(results{prc_10}{config_idx}(:));
prc_90_data = squeeze(results{prc_90}{config_idx}(:));

vals = [x_data, y_data, prc_10_data, prc_90_data];
good_idx = all(~isnan(vals),2);

xg = x_data(good_idx);
yg = y_data(good_idx);
prc10g = prc_10_data(good_idx);
prc90g = prc_90_data(good_idx);

prc_10_length = abs(yg - prc10g);
prc_90_length = abs(prc90g - yg);

[sorted_x, idx] = sort(xg);
upper = yg(idx) + prc_90_length(idx);
lower = yg(idx) - prc_10_length(idx);

patch_color = [1 0 0]; % Use your color
patch_alpha = 0.3;

prc_10_length = abs(y_data-prc_10_data);
prc_90_length = abs(y_data-prc_90_data);

hold on
fill([sorted_x; flipud(sorted_x)], [upper; flipud(lower)], [0,0,0], 'FaceAlpha', 0.3, 'EdgeColor', 'none')

errorbar(x_data,y_data,prc_10_length,prc_90_length, ...
    plotting_marker{config_idx}, ...
    'MarkerFaceColor',plotting_color(config_idx,:), ...
    'MarkerEdgeColor', 'k', ...
    'DisplayName',legend_label, ...
    LineStyle='-', ...
    Color='k', ...
    MarkerSize=8)

hold off
xlim([0 max(squeeze(results{x_num}{config_idx}(:)))+1])

if plot_legend == 1
    legend
end

box on

end