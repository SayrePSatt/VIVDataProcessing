function plot_psd_fn(results,x_num,y_num,z_num,config_idx,plot_legend,plotting_color,plotting_marker)
%This function takes in the results of f
% requency PSD contours for different
%tests and plots it on a spectrogram like plot
% results = results_ave;
% y_num = PSD_freq_norm_ave;
% z_num = PSD_norm_ave;
% x_num = 1;
springs_idx = 2;
% config_idx = 1;
num_datasets = length(y_num);
all_freqs = [];
for k = 1:num_datasets
    all_freqs = [all_freqs; y_num{k}(:)];
end
common_freq = linspace(min(all_freqs), max(all_freqs), 1000);

combined_data = NaN(length(common_freq), num_datasets);
for k = 1:num_datasets
    combined_data(:,k) = interp1(y_num{k}, z_num{k}, common_freq, 'linear', NaN);
end
x_data = cell2mat(squeeze(results{x_num}{config_idx}(:)));

imagesc(x_data, common_freq, combined_data)
set(gca, 'YDir', 'normal')
colorbar

xlabel('$U^*$')
colormap(flipud(gray))

clim([-3 0])
set(gca,'Colorscale','linear')
xlim([0 max(x_data)]);
ylim([0 6]);
% %%
% 
% x_data = cell2mat(squeeze(results{x_num}{config_idx}(:)));
% % [reps,~] = size(y_num);
% % x_data = repmat(x_data,reps);
% y_data = y_num(:,1);
% z_data = z_num;
% 
% x = x_data;
% y = y_data;
% z = z_data;
% 
% % Define regular grid to interpolate onto
% x_grid = linspace(min(x), max(x), 100);   % Adjust points for resolution
% y_grid = linspace(min(y), max(y), 1000);
% [X, Y] = meshgrid(x_grid, y_grid);
% 
% % Interpolate z onto the grid
% Z_grid = griddata(x, y, z, X, Y, 'linear');
% 
% imagesc(x_grid, y_grid, Z_grid)
% set(gca, 'YDir', 'normal')    % So y increases upwards
% colorbar
% % pcolor(x_data,y_data,z_data_mask);
% % shading interp;
% % if plot_legend == 1
% %     colorbar
% % end
% xlabel('$U^*$')
% ylabel('$f^*$')
% colormap(flipud(gray))
% 
% clim([-3 0])
% set(gca,'Colorscale','linear')
% xlim([0 max(x_data)]);
% ylim([0 7]);
% % 
% % h = get(gca,'children');
% % uistack(h(7),'top')
end