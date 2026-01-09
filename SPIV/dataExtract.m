tic
clear;clc;

fprintf('Loading files...\n')
data_files = dir ('000D_1k_06.5/*.dat');

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

path = strcat(data_files(1).folder,'\',data_files(1).name);
data = readmatrix(path);

I = find(isnan(data(:,1)));

% this loop is importatnt as some versions to MATLAB has a NaN at the top
% column
if ~isempty(I)
    x = data(I(end)+1:end,1)-trailingEdgeLocation;
    y = data(I(end)+1:end,2);
else
    x = data(:,1)-trailingEdgeLocation;
    y = data(:,2);
end

fprintf('Organizing files...\n')
u= zeros(length(x),length(data_files));
v= zeros(length(x),length(data_files));
w= zeros(length(x),length(data_files));
isValid = zeros(length(x),length(data_files));
parfor i = 1:length(data_files)
    path = strcat(data_files(i).folder,'\',data_files(i).name);
    data = readmatrix(path);
    if ~isempty(I)
        u(:,i) = data(I(end)+1:end,3);
        v(:,i) = data(I(end)+1:end,4);
        w(:,i) = data(I(end)+1:end,5);
        isValid(:,i) = data(I(end)+1:end,6);
    else
        u(:,i) = data(:,3);
        v(:,i) = data(:,4);
        w(:,i) = data(:,5);
        isValid(:,i) = data(:,6);
    end
end

fprintf('Saving data...\n')
save('vectorFields000D_1k_06.5.mat','x','y','u','v','w','isValid')
toc
