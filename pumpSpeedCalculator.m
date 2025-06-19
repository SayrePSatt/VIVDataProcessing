function [velocityUser] = pumpSpeedCalculator(requiredPumpSpeed)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% tic

% fprintf('Loading files...\n')
data_files = dir ('data/*.dat');

cwd = pwd;
d = dir(cwd);
% remove all files (isdir property is 0)
subFolder = d([d(:).isdir]) ;
% remove '.' and '..'
subFolder = subFolder(~ismember({subFolder(:).name},{'.','..'}));

trailingEdgeLocation = 0;  %in mm
chord = 90;    %in mm
% xLocation = [0.05 0.1 0.25 0.4 0.5 0.75 1]*chord;
% xLocationString = string(xLocation/chord);
% rho = 1000;
% V_TE_max = [0.222 0.222 0.222 0.222 0.222 0.0411 0.0411 0.0411];

% path = strcat(data_files(1).folder,'\',data_files(1).name);
% data = readmatrix(path);

frequency = [2 3 5 8 10 12 15];

% fprintf('Organizing files...\n')
% u= zeros(length(x),length(data_files));
% v= zeros(length(x),length(data_files));
for i = 1:length(data_files)
    path = strcat(data_files(i).folder,'\',data_files(i).name);
    data = readmatrix(path);
    I = find(isnan(data(:,1)));

    if ~isempty(I)
        x = data(I(end)+1:end,1);
        y = data(I(end)+1:end,2);
        u(:) = data(I(end)+1:end,3);
        v(:) = data(I(end)+1:end,4);
    else
        x = data(:,1);
        y = data(:,2);
        u(:) = data(:,3);
        v(:) = data(:,4);
    end

    VelocityTemp = u(u~=0);
    Velocity(i) = mean(VelocityTemp);
    VelocityStd(i) = std(VelocityTemp);
    clear x y u v VelocityTemp
end

% fprintf('Properties of the linear regression model fitted to the PIV data:');
mdl = fitlm(frequency,Velocity);
% fprintf(sprintf('\n\nPump Frequency of the channel flow for pump frequency of %d (m/s):',requiredPumpSpeed));
velocityUser = predict(mdl,requiredPumpSpeed);

% figure
% errorbar(frequency,Velocity,VelocityStd)
% toc

end