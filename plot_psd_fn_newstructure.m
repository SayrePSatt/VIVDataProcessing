function plot_psd_fn(results,x_num,y_num,z_num,config_idx,plot_legend,plotting_color,plotting_marker)
%This function takes in the results of f
% requency PSD contours for different
%tests and plots it on a spectrogram like plot

x_data = squeeze(results{x_num}{config_idx}(:));
y_data = y_num;
z_data = z_num;
z_data(z_data == 0) = NaN;

lower_limit = 10;
upper_limit = 60;

z_data_mask = z_data;
% z_data_mask(z_data_mask<lower_limit) = NaN;

pcolor(x_data,y_data,z_data_mask);
shading interp;
if plot_legend == 1
    colorbar
end
xlabel('$U^*$')
ylabel('$f^*$')
colormap(flipud(gray))

clim([-5 0])
set(gca,'Colorscale','linear')
xlim([0 max(x_data)]);
ylim([0 2]);
% 
% h = get(gca,'children');
% uistack(h(7),'top')
end