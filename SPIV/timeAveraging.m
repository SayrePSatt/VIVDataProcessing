clear;clc;close all;

tic 
R = [linspace(0,1,100) ones(1,43) linspace(1,1,100)];
G = [linspace(0,1,100) ones(1,43) linspace(1,0,100)];
B = [linspace(1,1,100) ones(1,43) linspace(1,0,100)];
cmap = [R(:),G(:),B(:)];  %// create colormap


chord = 0.0889;
VTEMax = 0.13;

% xLocation = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.6 1.7 1.8 1.9 2]*chord;
% xLocation = [0.05 0.25 0.4 0.75 1 1.5 2]*chord;
%% PIV data analysis

each_field=["vectorFields000D_1k_14.5.mat" "vectorFields040D_1k_14.5.mat"];

f = figure('units','inch','position',[1,1,8,6]);
tiledlayout(length(each_field),1,'TileSpacing', 'compact', 'Padding', 'none')

for each_field=["vectorFields000D_1k_14.5.mat" "vectorFields040D_1k_14.5.mat"]
nexttile
box on
xticklabels([]);
load(each_field)
x=(x/1000);
y=(y/1000);
[X,Y] = uniqueXYExtraction(x,y);
[XX,YY] = meshgrid(X,Y);

uAvg = mean(u,2);
vAvg = mean(v,2);
vFluc = mean(v-vAvg,2);

UmagAvg = sqrt(uAvg.^2+vAvg.^2);


ux = (reshape(uAvg,length(X),length(Y)))';
uy = (reshape(vAvg,length(X),length(Y)))';


% uyFluc = (reshape(vFluc,length(X),length(Y)))';
Umag = (reshape(UmagAvg,length(X),length(Y)))';
vorticity = curl(XX,YY,ux,uy);
vorticityStar = vorticity*chord/VTEMax;

% sigma = 3; % Standard deviation for the Gaussian filter
% h_gaussian = fspecial('gaussian', [25 25], sigma);
% vorticityStar = imfilter(vorticityStar, h_gaussian, 'replicate');

% [startX,startY] = meshgrid(-0.1:0.01:0.05,-0.06:0.01:0.06);

skip = 4;
% f = figure('units','inch','position',[1,1,8,6]);
% ax = gca;
% ax.LineWidth = 4;
% ax.FontSize = 18;
% ax.FontWeight = 'bold';
hold on
% contourf(X,Y,Gamma1)
contourf(X/chord,Y/chord,vorticityStar,500,'LineColor','none','edgecolor','none','LevelList',-11:0.005:11)
clim([-2 2])
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
xlim([-1 1])
ylim([-0.5 0.5])
ax1 = gca;
ax1.XMinorTick = 'on';
ax1.XAxis.MinorTickValues = -1:0.1:1;
ax1.XAxis.TickValues = [-1 0 1];
ax1.YMinorTick = 'on';
ax1.YAxis.MinorTickValues = -0.5:0.1:0.5;
ax1.YAxis.TickValues = [-0.5 0 0.5];
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
center_x = 0.5;
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

end
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

%% Graphics export
exportgraphics(f,['SPIV_figures\' 'vectorFields000D_1k_14.5.pdf'],'Resolution',600);
exportgraphics(f,['SPIV_figures\' 'vectorFields000D_1k_14.5.png'],'Resolution',300);
