clear all
close all
clc

load("pumpFit_freq2velo.mat");

topfolder = 'D:\EFDL\NEW_VIVmaster\tandemSphereData\velo_fields\02_24_2026\outoffoler\';
[temp, processing_file_dir] = uigetfile(topfolder,'*.vc7','MultiSelect','on');

if class(temp) == 'char'
    processing_files = cellstr(temp);
    num_figs = 1;
else
    processing_files = temp;
    num_figs = numel(processing_files);
end

tic 
% R = [linspace(0,1,100) ones(1,43) linspace(1,1,100)];
% G = [linspace(0,1,100) ones(1,43) linspace(1,0,100)];
% B = [linspace(1,1,100) ones(1,43) linspace(1,0,100)];
% cmap = [R(:),G(:),B(:)];  %// create colormap

reds = [251 226 211;
        246 178 148;
        224 109 84]/255;       % darkest red

whites = [1 1 1];     % white

blues = [219 234 242;   % lightest blue
         156 199 223;
         33 102 172]/255;      % darkest blue

% colors = [224 109 84;
%           246 178 148;
%           251 226 211;   % darkest red
%           255 255 255;    % white
%           219 234 242;   % lightest blue
%           156 199 223;
%           33 102 172]/255;      % darkest blue
colors = [
    165 32 38;      % new darkest red
    165 32 38;      % new darkest red
    194 70 61;      % deep red (added)
    194 70 61;      % deep red (added)
    224 109 84;     % mid red
    224 109 84;     % mid red
    246 178 148;    % lighter red
    246 178 148;    % lighter red
    251 226 211;    % lightest red
    251 226 211;    % lightest red
    255 255 255;    % white
    219 234 242;    % lightest blue
    219 234 242;    % lightest blue
    156 199 223;    % lighter blue
    156 199 223;    % lighter blue
    89 150 197;     % mid blue (added)
    89 150 197;     % mid blue (added)
    33 102 172;     % darkest blue
    33 102 172;     % darkest blue
    12 44 132       % new deeper blue
    12 44 132       % new deeper blue
] / 255;

% colors = turbo(11);
% colors(5,:) = [1 1 1];

colors = flipud(colors);
cmap = colors;

%%
make_subplot = 1;
if make_subplot == 1
    subplot_fig = figure('Position',[188 194 1519 721]);
    for jj = 1:numel(processing_files)
        dists(jj) = str2double(processing_files{jj}(1:3));
        subbin(jj) = str2double(processing_files{jj}(45:46));
    end
    num_col = numel(unique(dists));
    num_row = numel(unique(subbin));
    t_fig = tiledlayout(num_row,num_col);
    t_fig.TileSpacing = 'none';
    t_fig.Padding = 'none';
    t_fig.TileIndexing = 'columnmajor';
end


plot_quiver = 0;
skip = 10;
scale = 1;
max_vorticity = 1.5;
smoothing_window = 7; %7 for vorticity
sigma = 1.1;
N=2;
value_to_plot = 1; %1 for vorticity, 2 for u, 3 for v, 4 for w
sliding_ave = ones(smoothing_window,smoothing_window)/smoothing_window^2;
min_size = 1.2/100; %Minimum window size for any vortex to be saved

xlim_vals = [-1.5 1.5];
ylim_vals = [-0.7 0.7];

f = figure;
for ii=1:length(processing_files)
% pause
% f = figure('units','inch','position',[1,1,8,6]);
close(f(isvalid(f)))
diameter = str2double(processing_files{ii}(13:14))/1000;
f_pump = str2double(processing_files{ii}(24:28));
[U U_68_temp] = predict(mdl,f_pump,Alpha=0.05);
distance = str2double(processing_files{ii}(1:3))/10;
u_star = (processing_files{ii}(36:39));
f = figure('units','inch','position',[1,1,8,6]);
bin = str2double(processing_files{ii}(45:46));
nBins = str2double(processing_files{ii}(53:54));
% tiledlayout(length(each_field),1,'TileSpacing', 'compact', 'Padding', 'none')

%%
rho = 998;
C_A = 0.5;     %Added mass coefficient
St = 0.19;
St_68 = 0.005;

diagnose = false;

markers = ['s' 'd' '*'];

field = readimx(fullfile(processing_file_dir,processing_files{ii}));
fieldFrame = field.Frames{1};
%%
D = create2DVec(fieldFrame);
D.W = fieldFrame.Components{5}.Planes{1};

X = D.X(:,1)/1000;
Y = D.Y(1,:)/1000;
[XX, YY] = meshgrid(X,Y);

Nx = N*size(X,1);
Ny = N*size(Y,2);
[XXq, YYq] = meshgrid(linspace(min(X),max(X),Nx), linspace(min(Y),max(Y),Ny)); % N larger than original size

if value_to_plot == 1    

    vorticity = curl(X,Y,D.U',D.V');
    vorticityStar = vorticity*diameter/U;
    
    
    plot_parameter = vorticityStar;
    gaussian_kernel = fspecial('gaussian', smoothing_window, sigma);
    plot_parameter = conv2(vorticityStar,sliding_ave,'same');

    max(abs(plot_parameter),[],"all")
    % 
    
    % alpha = 0.1;           % Filter strength; try between 0.01 and 0.2 for tuning
    % nIterations = 3;      % How many times to apply filter? (e.g., 5-50)
    % 
    % for k = 1:nIterations
    %     % Compute discrete Laplacian
    %     laplacian = -4*u + ...
    %                 circshift(u, [1, 0]) + ... % up
    %                 circshift(u, [-1, 0]) + ... % down
    %                 circshift(u, [0, 1]) + ... % right
    %                 circshift(u, [0, -1]);     % left
    % 
    %     % Update with Laplacian (diffusion smoothing)
    %     u = u + alpha * laplacian;
    % end
    % % sigma = 3; % Standard deviation for the Gaussian filter
    % % h_gaussian = fspecial('gaussian', [25 25], sigma);
    % % u = imfilter(vorticityStar, h_gaussian, 'replicate');
    % 
    % % [startX,startY] = meshgrid(-0.1:0.01:0.05,-0.06:0.01:0.06);
    
    % skip = 4;
    % f = figure('units','inch','position',[1,1,8,6]);
    % ax = gca;
    % ax.LineWidth = 4;
    % ax.FontSize = 18;
    % ax.FontWeight = 'bold';
    hold on

    plot_para_smooth = interp2(XX, YY, plot_parameter, XXq, YYq, 'cubic');
    plot_para_smooth(abs(plot_para_smooth)<0.2) = 0;
    
    % levels = [-3 -2.1429 -1.2677 -1.2677 2.1429 3];
    plot_para_smooth(plot_para_smooth>=max_vorticity) = max_vorticity;
    plot_para_smooth(plot_para_smooth<=-max_vorticity) = -max_vorticity;
    BW = abs(plot_para_smooth) > 0.2;
    CC = bwconncomp(BW); 
    region_area = cellfun(@numel, CC.PixelIdxList);  % Number of pixels (area) for each region
    image_area = prod(CC.ImageSize);

    % Sort regions by area descending
    [~, sort_idx] = sort(region_area, 'descend');

    % Create mask for the two largest regions
    mask = zeros(size(BW));  % Initialize mask
    numRegionsToKeep = min(4, length(sort_idx));  % In case fewer than 2 regions exist

    for k = 1:numRegionsToKeep
        if numel(CC.PixelIdxList{sort_idx(k)})>min_size*image_area
            mask(CC.PixelIdxList{sort_idx(k)}) = 1;  % Set pixels for top regions to 1
        end
    end

    %Keep only the two largest vortices field
    plot_para_smooth = plot_para_smooth .* mask;
    
    cbar_name = '$\omega^*$';
    plot_name = '_vorticity';
elseif value_to_plot == 2
    plot_parameter = D.U'/U;
    cbar_name = '$u/U_{\infty}$';
    plot_name = '_u';
    plot_parameter = conv2(plot_parameter,sliding_ave,'same');
    plot_para_smooth = interp2(XX, YY, plot_parameter, XXq, YYq, 'cubic');
elseif value_to_plot == 3
    plot_parameter = D.V'/U;
    cbar_name = '$v/U_{\infty}$';
    plot_name = '_v';
    plot_parameter = conv2(plot_parameter,sliding_ave,'same');
    plot_para_smooth = interp2(XX, YY, plot_parameter, XXq, YYq, 'cubic');
elseif value_to_plot == 4
    plot_parameter = D.W'*2/(100*U); %/100 because for some reason thru plane velocity reported in cm/s
    cbar_name = '$w/U_{\infty}$';
    plot_name = '_w';
    plot_parameter = conv2(plot_parameter,sliding_ave,'same');
    plot_para_smooth = interp2(XX, YY, plot_parameter, XXq, YYq, 'cubic');
end

U_interp = interp2(XX, YY, D.U', XXq, YYq, 'cubic');
V_interp = interp2(XX, YY, D.V', XXq, YYq, 'cubic');

levels = [-max_vorticity:2*max_vorticity/7:max_vorticity];
levels(4:5) = [];
levels = linspace(-max_vorticity,max_vorticity,length(colors));

% contourf(X/diameter,Y/diameter,u,levels,'LineColor','none','edgecolor','k')%,'LevelList',-max_vorticity:max_vorticity/7:max_vorticity)
figure(f)
contourf(XXq/diameter,YYq/diameter,plot_para_smooth,levels,'LineColor','none')%,'LevelList',-max_vorticity:max_vorticity/7:max_vorticity)
clim([-max_vorticity max_vorticity])
colormap(cmap)

if plot_quiver == 1
    quiver(XXq(1:skip:end,1:skip:end)/diameter,YYq(1:skip:end,1:skip:end)/diameter,U_interp(1:skip:end,1:skip:end)/U,V_interp(1:skip:end,1:skip:end)/U,'k')
end
% quiver(X(1:skip:end)/chord,Y(1:skip:end)/chord,ux(1:skip:end,1:skip:end),uy(1:skip:end,1:skip:end),1,'k','LineWidth',1)
% h = stream2(X/chord,Y/chord,ux,uy,startX,startY);
% lineobj = streamline(h);
% set(lineobj, 'Color', [0.5 0.5 0.5])
% scatter(xLoc,yLoc,60,'filled',"yellow")
% scatter(mean_x,mean_y,40,'filled',"green")
hold off

% set(gca,'FontWeight','Bold')

%     grid on
ylabel('$z/D$','Rotation',0)
% ,'Position',[-120 170]
% xlabel('t (s)','FontAngle','italic')
% xlabel('y/D')
axis equal
xlim(xlim_vals)
ylim(ylim_vals)
ax1 = gca;
ax1.XMinorTick = 'on';
ax1.XAxis.MinorTickValues = -1.5:0.1:1.5;
ax1.XAxis.TickValues = [-1 0 1];
ax1.YMinorTick = 'on';
ax1.YAxis.MinorTickValues = -0.7:0.1:0.7;
ax1.YAxis.TickValues = [-0.7 0 0.7];
% xticks([-1:0.1:1])
% set(ax1,'Xticklabel',[-1 0 1])
set(ax1,'FontWeight','Bold','fontsize',20,'FontName','Times','LineWidth',2,'TickLength',[0.02 0.02],'Layer','top')
% saveas(f,'Re3700SmoothTimeAveraged.png');

% save 'smoothFlowField2Re3700.mat' 'X' 'Y' 'ux' 'uy' 'vorticityStar'

% save 'smoothFlowField3Re3700.mat' 'X' 'Y' 'ux' 'uy' 'vorticityStar'

% [yPoints,uSlice] = extractSliceVelocities(XX,YY,Y,ux,xLocation);
% [~,vSlice] = extractSliceVelocities(XX,YY,Y,uy,xLocation);
% [~,vSliceFluc] = extractSliceVelocities(XX,YY,Y,uyFluc,xLocation);
% [~,UmagSlice] = extractSliceVelocities(XX,YY,Y,Umag,xLocation);
% [meanUx,maxUx,maxLoc,maxUy] = meanNmaxU(xLocation,yPoints,uSlice,vSlice);
% [meanUmag,maxUmag,~,~] = meanNmaxU(xLocation,yPoints,UmagSlice,vSlice);
% save 'smoothVelocityExtractRe3700Updated.mat' 'yPoints' 'uSlice' 'meanUx' 'maxUx' 'maxLoc' 'UmagSlice' 'meanUmag' 'maxUmag'
% save 'smoothThrustRe20000.mat' 'yPoints' 'uSlice' 'vSlice' 'maxUx' 'maxLoc' 'maxUy' 'vSliceFluc'
% save 'timeAveragedVortexTestRe10000Smooth.mat' 'X' 'Y' 'uAvg' 'vAvg'
% save 'timeAveragedThrustTestRe10000Smooth.mat' 'X' 'Y' 'uAvg' 'vAvg' 'vFluc'

toc

%% Circle
% Your desired circle parameters
hold on
u_red = str2double(u_star);
if distance == 0
    if u_red == 8.5
        scale = 0.7;
    elseif u_red == 14.5
        scale = 0.6;
    end
elseif distance == 1.5
    if u_red == 8.5
        scale = 0.6;
    elseif u_red == 14.5
        scale = 1.1;
    end
elseif distance == 2.0
    if u_red == 8.5
        scale = 0.5;
    elseif u_red == 14.5
        scale = 1.05;
    end
elseif distance == 2.5
    u_red = 0.6;
elseif distance == 4.0
    if u_red == 8.5
        scale = 0.7;
    elseif u_red == 14.5
        scale = 0.95;
    end
end

pos_x_sphere = scale*round(sin(2*pi*(bin-1)/nBins),2);
velo_x_sphere = round(cos(2*pi*(bin-1)/nBins),2);
center_x = 1*pos_x_sphere;
center_y = 0;
diameter = 1;

% Calculate the position vector [x_left, y_bottom, width, height]
if velo_x_sphere<0
    x_arrow_pos = center_x + diameter/4;
    arrow_dir = -1*diameter;
elseif velo_x_sphere>0
    x_arrow_pos = center_x - diameter/4;
    arrow_dir = 1*diameter;
elseif velo_x_sphere==0
    x_arrow_pos = center_x;
    arrow_dir = 0;
end
pos = [center_x - diameter/2, center_y - diameter/2, diameter, diameter];

% Draw the dotted circle
rectangle('Position', pos, 'Curvature', [1, 1], ...
          'LineStyle', '--',               ... % Dotted line
          'EdgeColor', [128 128 128]/255,              ... % Color: black
          'LineWidth', 2);                   % Optional: make it visually bold
arrow_length = arrow_dir*0.4;  % Length of arrow (adjust as needed)
quiver(x_arrow_pos, center_y, arrow_length, 0, 0,'Color', [128 128 128]/255, 'LineWidth', 2, 'MaxHeadSize', 2);

hold off;


%% Figure labeling
figure(f)
hold off
set(gca,'TickLabelInterpreter','latex')
colormap(cmap)
% daspect([8 8 1])
h = colorbar;
drawnow
h.TickLabelInterpreter = 'latex';
h.Label.Interpreter = 'latex';
h.Label.String = cbar_name;
h.Label.Rotation = 0;
h.Location = 'eastoutside';
h.FontWeight = 'normal';

% h.Position = [0.88 0.3000 0.0208 0.5000];
xlabel('$y/D$')
xticklabels({'-1', '0', '1'})

if distance == 0
    dist_disp = ', Isolated';
else
    dist_disp = [', $L^*=$' num2str(distance)];
end

title(['$U^*=$' num2str(str2double(u_star)) dist_disp])
box on


if make_subplot == 1
    figure(subplot_fig)
    nexttile
    % subplot(num_row,num_col,ii)
    axis equal
    box on
    copyobj(get(findobj(f,'Type','Axes'),'Children'),gca)
    set(gca,'TickLabelInterpreter','latex')
    xlim(xlim_vals)
    ylim(ylim_vals)
    clim([-max_vorticity max_vorticity])
    colormap(cmap)
    % daspect([8 8 1])
    % h_sp = colorbar;
    if ii == num_col*num_row
        h_sp = colorbar;
        h_sp.TickLabelInterpreter = 'latex';
        h_sp.Label.Interpreter = 'latex';
        h_sp.Label.String = cbar_name;
        h_sp.Label.Rotation = 0;
        h_sp.Layout.Tile = 'east';
        h_sp.FontWeight = 'normal';
        drawnow
    end
    
    if mod(ii,num_row) == 1
        title(['$U^*=$' num2str(str2double(u_star)) dist_disp])
    end
    
    ax1 = gca;
    ax1.XMinorTick = 'on';
    ax1.XAxis.MinorTickValues = xlim_vals(1):0.1:xlim_vals(2);
    ax1.XAxis.TickValues = [-1 0 1];
    if not(mod(ii,num_row) == 0)
        ax1.XAxis.TickLabels = {};
    else
        xlabel('$y/D$','Rotation',0)
    end
    
    ax1.YMinorTick = 'on';
    ax1.YAxis.MinorTickValues = ylim_vals(1):0.1:ylim_vals(2);
    ax1.YAxis.TickValues = [-0.5 0 0.5];
    if not(ii<=num_row)
        ax1.YAxis.TickLabels = {};
    else
        ylabel('$z/D$','Rotation',0)
    end
    % xticks([-1:0.1:1])
    % set(ax1,'Xticklabel',[-1 0 1])
    set(ax1,'FontWeight','Bold','fontsize',20,'FontName','Times','LineWidth',2,'TickLength',[0.02 0.02],'Layer','top')
    
    cur_row = mod(ii,num_row);
    if cur_row == 0
        cur_row = num_row;
    end
    
    cur_col = ceil(ii/num_row);
    subfig_label = [char(cur_col+96) '.' num2roman(cur_row) ')'];
    
    text(ax1, 0.04, 0.96, subfig_label, ...
        'Units', 'normalized', ...
        'FontWeight', 'bold', ...
        'FontSize', 20, ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left');
end

exportgraphics(f,fullfile([processing_file_dir '\figures\'], [processing_files{ii}(1:end-4) plot_name '.pdf']),'Resolution',600);
exportgraphics(f,fullfile([processing_file_dir '\figures\'], [processing_files{ii}(1:end-4) plot_name '.png']),'Resolution',300);

end
subplot_name = [processing_files{ii}(10:39) '_smooth_' num2str(smoothing_window) plot_name];
exportgraphics(subplot_fig,fullfile([processing_file_dir '\figures\'], [subplot_name '_subplot.png']),'Resolution',600);