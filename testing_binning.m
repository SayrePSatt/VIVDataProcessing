clear all
close all
clc
% Testing binning

t = 0:0.001:5;
omega = 2*pi*0.25;
disp = sin(omega*t);
velo = omega*cos(omega*t);

disp_norm = disp/max(disp);
velo_norm = velo/max(velo);
dispvelo_mag_norm = sqrt(disp_norm.^2+velo_norm.^2);
angle_norm = atan2(disp_norm,velo_norm);
angle_norm = wrapTo2Pi(angle_norm);

nBins = 24;
edges = linspace(0,2*pi,(nBins*2)+1);
edges = sort(edges); % make sure they are sorted
binned_data = discretize(angle_norm,edges);

jj = 2;
for ii = 2:2:(nBins*2)-1
    binned_data(binned_data==ii) = jj;
    binned_data(binned_data==ii+1) = jj;
    jj = jj+1
end
binned_data(binned_data==nBins*2) = 1;

unique(binned_data)

colors = lines(nBins);
color_matrix = colors(binned_data,:);
hold on
scatter(velo_norm,disp_norm,20,color_matrix,"filled")

plot(velo/max(velo),disp/max(disp))
axis equal
