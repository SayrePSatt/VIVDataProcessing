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


x_data = cell2mat(squeeze(results{x_num}{config_idx}(:)));
springs = cell2mat(squeeze(results{springs_idx}{config_idx}(:)));
uniq_springs = unique(springs);

x_data_1k = x_data(springs==uniq_springs(1))';
x_data_6k = x_data(springs == uniq_springs(2))';
% 
% x_data_1k = linspace(min(x_data_1k),max(x_data_1k),80);
% x_data_6k = linspace(min(x_data_6k),max(x_data_6k),15);

y_num_1k = y_num(springs == uniq_springs(1));
y_num_6k = y_num(springs == uniq_springs(2));
y_data_1k = y_num_1k{1}';
y_data_6k = y_num_6k{1}';

z_num_1k = z_num(springs == uniq_springs(1));
z_num_6k = z_num(springs == uniq_springs(2));
z_data_1k = cell2mat(z_num_1k);
z_data_6k = cell2mat(z_num_6k);

common_y = linspace(min([y_data_1k y_data_6k]), max([y_data_1k y_data_6k]), 1000);

Z1_interp = zeros(length(common_y), length(x_data_1k));
for col = 1:length(x_data_1k)
    Z1_interp(:,col) = interp1(y_data_1k, z_data_1k(:,col), common_y, 'linear', NaN);
end

Z2_interp = zeros(length(common_y), length(x_data_6k));
for col = 1:length(x_data_6k)
    Z2_interp(:,col) = interp1(y_data_6k, z_data_6k(:,col), common_y, 'linear', NaN);
end

x_combined = [x_data_6k x_data_1k];
Z_combined = [Z2_interp Z1_interp];

imagesc(x_combined, common_y, Z_combined)
set(gca, 'YDir', 'normal')
colorbar

xlabel('$U^*$')
ylabel('$f^*$')
colormap(flipud(gray))

clim([-3 0])
set(gca,'Colorscale','linear')
xlim([0 max(x_data)]);
ylim([0 7]);
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