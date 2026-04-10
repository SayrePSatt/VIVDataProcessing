%This code takes in a data file and bins the file, exporting the indicies
%of images to export from DaVIS
clear all
close all
clc

load("pumpFit_freq2velo.mat");
%% File Import
[processing_files, processing_file_dir] = uigetfile('*.dat','MultiSelect','on');
SPIV_folder = uigetdir(processing_file_dir);

exclude_zeros = '00.00';
zero_mask = ~contains(processing_files,exclude_zeros);
processing_files = processing_files(zero_mask);

for ii = 1:length(processing_files)
file = processing_files{ii};
file_dir = [processing_file_dir, file];
%% File Read
test_spec = regexp(file,'\d\d\dD_\d\dD');

metadata = table2array(readtable(file_dir,'Range','A12:F13'));
data = table2array(readtable(file_dir,'NumHeaderLines',14)); %Imports one file with corresponding data
% if data(end,1) > 225
%     data = data(100000:end,:);
% end

%% Extracting metadata and run specifications
rho = 998;
temp_loc = regexp(file,'\d\d.\d\d.dat');
f_pump = str2num(file(temp_loc:temp_loc+4));
[U U_68_temp] = predict(mdl,f_pump,Alpha=0.05);%/(1.117645);
U_68 = U-U_68_temp(1);
d_sph = metadata(:,1);
m_d = (4/3)*pi*(d_sph(1)/2)^3*rho+rho*0.005^2*pi*d_sph(1)/4; %Displaced mass
m = metadata(:,2);
f_nw = metadata(:,3);
f_na = metadata(:,5);
zeta = metadata(:,6);

m_d = (4/3)*pi*(d_sph(1)/2)^3*rho+rho*0.005^2*pi*d_sph(1)/4;
m_star = m(1)/m_d;
C_A = ((f_na(1)./f_nw(1)).^2-1)*m_star;
m_a = C_A*m_d;
omega_na = 2*pi*f_na(1);
k = m*omega_na.^2; %5.375; %(f_n(1)*2*pi)^2*m
c = 2*m*omega_na*zeta(1);
mass_damp = (m_star+C_A)*zeta(1);

time = data(:,1);
f_s = 1/(time(2)-time(1));
dt = 1/f_s;
encoder = data(:,2);

%% Filtering
% subplot(3,1,1)
n = length(time);
fhat = fft(encoder, n); % Compute the fast fourier transform
PSD = fhat.*conj(fhat)/n; % Power Spectrum(power per freq)

freq = 1/(dt*n)*(0:n); % Create x-qxis of frequencies in Hz
L = 1:floor(n/2); %only plot the first half of freqs
f_c = 2;   %Cutoff frequency
[b,a] = butter(4,f_c/(f_s/2),"low");

encoder_filt = filtfilt(b,a,encoder);

%% Data Processing
velo = FivePointDiff(encoder_filt,f_s)';
% rms(velo)
% velo = filtfilt(b,a,velo);
acc = FivePointDiff(velo,f_s)';
data_length = length(acc);
time = time(1:data_length);
encoder_filt = encoder_filt(1:data_length);
velo =  velo(1:data_length);

%% Zero correction from encoder drift
[peak_velo zero_pos] = findpeaks(velo);
zero_pos = zero_pos(1:end-1);
zero_fit = polyfit(time(zero_pos),encoder_filt(zero_pos),1);
% nodeIndices = find(acc(1:end-1).*acc(2:end) < 0);
% nodeIndices = nodeIndices(1:end-1);
% 
% figure
% plot(time(zero_pos),encoder_filt(zero_pos))
% hold on
% plot(time(zero_pos),polyval(zero_fit,time(zero_pos)))

encoder_filt_zeroed = encoder_filt - polyval(zero_fit,time);
%Maybe think about recalculating velocity. Calculations of acceleration
%shouldn't change
% 
% figure
% plot(time(zero_pos),encoder_filt_zeroed(zero_pos))
% hold on
% plot(time(nodeIndices),encoder_filt(nodeIndices))

%% Getting Image Capture Locations
% Detecting changes from 0 to 1
clear sync ttl start_idx
sync = data(1:data_length,4);
ttl = data(1:data_length,5);
start_idx = find(sync(1:end-1) == 0 & sync(2:end) == 1) + 1;
for ii = 1:length(start_idx)
    jj = 1;
    while(ttl(start_idx(ii))>ttl(start_idx(ii)+jj)-1) %-1 because of weird timing signal
        jj = jj+1;
    end
    start_idx(ii,2) = start_idx(ii,1)+jj;
end

start_idx = start_idx(end,:); %takes only the last sync signal
clear end_idx
end_idx = find(sync(1:end-1) == 1 & sync(2:end) == 0);
for ii = 1:length(end_idx)
    jj = 1;
    while(ttl(end_idx(ii))==ttl(end_idx(ii)-jj)) %-1 because of weird timing signal
        jj = jj+1;
    end
    end_idx(ii,2) = end_idx(ii,1)-jj;
end

end_idx = end_idx(end,:);


%%Below is a code for changing the sync occurance if sync occured after the
%%trigger
% for ii = 1:length(change_idx)
%     if data(change_idx(ii),5) == data(change_idx(ii)-1,5)
%         change_idx(ii)
%     end
%     while data(change_idx(ii),5) == data(change_idx(ii)-1,5)
%         disp("changed");
%         change_idx(ii) = change_idx(ii)+1;
%     end
% end
sync_disp_norm = encoder_filt_zeroed(start_idx)/max(encoder_filt_zeroed);
sync_velo_norm = velo(start_idx)/max(velo);
%% Plotting binning
disp_norm = encoder_filt_zeroed/max(encoder_filt_zeroed);
velo_norm = velo/max(velo);
dispvelo_mag_norm = sqrt(disp_norm.^2+velo_norm.^2);
used_idx = dispvelo_mag_norm >= 0; %abs(dispvelo_mag_norm) > 0.7;
angle_norm = atan2(disp_norm,velo_norm);
angle_norm = wrapTo2Pi(angle_norm);

nBins = 24;
bin_forprocessing = [20 2 8 14];

edges = linspace(0,2*pi,(nBins*2)+1);
edges = sort(edges); % make sure they are sorted
binned_data = discretize(angle_norm,edges);

jj = 2;
for ii = 2:2:(nBins*2)-1
    binned_data(binned_data==ii) = jj;
    binned_data(binned_data==ii+1) = jj;
    jj = jj+1;
end
binned_data(binned_data==nBins*2) = 1;

for kk = 1:length(bin_forprocessing)
    processing_idx = find(binned_data(start_idx(2):end_idx(2)) == bin_forprocessing(kk));
    ttl_synced = ttl - ttl(start_idx(2)) + 1;
    % ttl_synced(start_idx(2)+1:2:end_idx(2)+1) = ttl(start_idx(2)+1:2:end_idx(2)+1);
    ttl_mask = ttl_synced(start_idx(2):end_idx(2));
    ttl_processing = unique(ttl_mask(processing_idx));
    ttl_processing = ttl_processing/2;
    ttl_2ndremoved = mod(ttl_processing,1) == 0;
    ttl_processing = (ttl_processing(ttl_2ndremoved))';
    
    %% Saving Bin file
    foldername = [SPIV_folder '\bin' num2str(bin_forprocessing(kk)) '_numbins' num2str(nBins)]; % Your folder name
    
    % Make sure the folder exists
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    
    bin_filename = ['bin' num2str(bin_forprocessing(kk)) '_' file(1:end-3) 'txt'];
    
    filename = fullfile(foldername,bin_filename);
    
    writematrix(ttl_processing, filename, 'Delimiter', ';');
    % fileID = fopen(filename,'w');
    % fprintf(fileID,'%d;',ttl_processing);
    % fclose(fileID);
    % unique(binned_data);
    % 
    % colors = lines(nBins);
    % color_matrix = colors(binned_data(used_idx),:);
    % close all
    % hold on
    % data_s = scatter(velo_norm(used_idx),disp_norm(used_idx),10,color_matrix,"filled")
    % data_s.HandleVisibility = 'off';
    
    % plot(velo_norm,disp_norm)
    % s1 = scatter(sync_velo_norm(:,1),sync_disp_norm(:,1),"filled","o",'k')
    % s1.DisplayName = '1st Image Pair';
    % s2 = scatter(sync_velo_norm(:,2),sync_disp_norm(:,2),"filled","^",'k')
    % s2.DisplayName = '2nd Image Pair';
    % legend
    % 
    % xlabel('$\dot{y}/\dot{y}_{max}$')
    % ylabel('$y/y_{max}$')
    % axis equal
end
end 