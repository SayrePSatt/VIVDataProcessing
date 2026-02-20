clear all
close all
clc

load psd_freq_norm_ave.mat
load psd_norm_ave.mat
load results_ave.mat


results = results_ave;
y_num = PSD_freq_norm_ave;
z_num = PSD_norm_ave;
x_num = 1;
config_idx = 1;
multi_spring = 1;

x_data = cell2mat(squeeze(results{x_num}{config_idx}(:)));
springs = [6 6 6 6 1 1 1 1 1 1 1 1 1];
uniq_springs = unique(springs);

x_data_1k = x_data(springs==uniq_springs(1))';
x_data_6k = x_data(springs == uniq_springs(2))';

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


