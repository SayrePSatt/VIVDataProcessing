clear all
close all
clc

load("pumpFit_freq2velo.mat");

topfolder = 'E:\EFDL\NEW_VIVmaster\tandemSphereData\SPIV\SPIV_forprocessing\';
[processing_files, processing_file_dir] = uigetfile(topfolder,'*.vc7','MultiSelect','on');

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

colors = [224 109 84;
          246 178 148;
          251 226 211;   % darkest red
          255 255 255;    % white
          219 234 242;   % lightest blue
          156 199 223;
          33 102 172]/255;      % darkest blue

colors = flipud(colors);
cmap = colors;
%%

max_vorticity = 2;

for ii=1:length(processing_files)

diameter = str2double(processing_files{ii}(13:14))/1000;
f_pump = str2double(processing_files{ii}(18:22));
[U U_68_temp] = predict(mdl,f_pump,Alpha=0.05);
distance = str2double(processing_files{ii}(1:3))/10;
u_star = (processing_files{ii}(30:33));
f = figure('units','inch','position',[1,1,8,6]);
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
[XX YY] = meshgrid(X,Y);

vorticity = curl(X,Y,D.U',D.V');
vorticityStar = vorticity*diameter/U;


u = vorticityStar;
u(abs(u)<0.2) = 0;
alpha = 0.06;           % Filter strength; try between 0.01 and 0.2 for tuning
nIterations = 5;      % How many times to apply filter? (e.g., 5-50)

for k = 1:nIterations
    % Compute discrete Laplacian
    laplacian = -4*u + ...
                circshift(u, [1, 0]) + ... % up
                circshift(u, [-1, 0]) + ... % down
                circshift(u, [0, 1]) + ... % right
                circshift(u, [0, -1]);     % left
    
    % Update with Laplacian (diffusion smoothing)
    u = u + alpha * laplacian;
end
% sigma = 3; % Standard deviation for the Gaussian filter
% h_gaussian = fspecial('gaussian', [25 25], sigma);
% u = imfilter(vorticityStar, h_gaussian, 'replicate');

% [startX,startY] = meshgrid(-0.1:0.01:0.05,-0.06:0.01:0.06);
BW = abs(u) > 0.3;
CC = bwconncomp(BW); 
region_area = cellfun(@numel, CC.PixelIdxList);  % Number of pixels (area) for each region

% Sort regions by area descending
[~, sort_idx] = sort(region_area, 'descend');

% Create mask for the two largest regions
mask = zeros(size(BW));  % Initialize mask
numRegionsToKeep = min(2, length(sort_idx));  % In case fewer than 2 regions exist

for k = 1:numRegionsToKeep
    mask(CC.PixelIdxList{sort_idx(k)}) = 1;  % Set pixels for top regions to 1
end

% Keep only the two largest vortices in your field
u = u .* mask;

skip = 4;
% f = figure('units','inch','position',[1,1,8,6]);
% ax = gca;
% ax.LineWidth = 4;
% ax.FontSize = 18;
% ax.FontWeight = 'bold';
hold on
levels = [-max_vorticity:2*max_vorticity/7:max_vorticity];
levels(4:5) = [];
levels = linspace(-max_vorticity,max_vorticity,7);
% levels = [-3 -2.1429 -1.2677 -1.2677 2.1429 3];
u(u>=max_vorticity) = max_vorticity;
u(u<=-max_vorticity) = -max_vorticity;
contourf(X/diameter,Y/diameter,u,levels,'LineColor','none','edgecolor','k')%,'LevelList',-max_vorticity:max_vorticity/7:max_vorticity)
clim([-max_vorticity max_vorticity])
colormap(cmap)
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
xlim([-1.5 1.5])
ylim([-0.7 0.7])
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
center_x = 0;
center_y = 0;
diameter = 1;

% Calculate the position vector [x_left, y_bottom, width, height]
pos = [center_x - diameter/2, center_y - diameter/2, diameter, diameter];

% Draw the dotted circle
rectangle('Position', pos, 'Curvature', [1, 1], ...
          'LineStyle', '--',               ... % Dotted line
          'EdgeColor', 'k',              ... % Color: black
          'LineWidth', 3.5);                   % Optional: make it visually bold
arrow_length = diameter * 0.4;  % Length of arrow (adjust as needed)
quiver(center_x, center_y, arrow_length, 0, 0, 'k', 'LineWidth', 4, 'MaxHeadSize', 2);

hold off;


%% Figure labeling
figure(f)
colormap(cmap)
% daspect([8 8 1])
h = colorbar;
drawnow
h.Label.String = '\omega^*';
h.Label.Rotation = 0;
h.Location = 'eastoutside';
h.FontWeight = 'normal';
h.Position = [0.88 0.3000 0.0208 0.5000];
xlabel('$y/D$')
xticklabels({'-1', '0', '1'})

if distance == 0
    dist_disp = 'Single';
else
    dist_disp = num2str(distance);
end

title(['$U^*$=' num2str(u_star) ' $L^*=$' dist_disp])

exportgraphics(f,fullfile(processing_file_dir, [processing_files{ii} '_vorticity.pdf']),'Resolution',600);
exportgraphics(f,fullfile(processing_file_dir, [processing_files{ii} '_vorticity.png']),'Resolution',300);

end