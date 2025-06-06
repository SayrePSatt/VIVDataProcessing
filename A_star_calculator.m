clear all
close all
clc

m = 10;
for i=1:3
    zero = csvread("test"+i+"/testVIV_"+i+"_00.00.csv");
    disp_zero_time = zero(:,1);
    disp_zero_pos = zero(:,2);
    disp_zero_mean = mean(disp_zero_pos);

pump_freq = ["03.00" "03.75" "04.50" "05.25" "06.00" "06.75" "07.50" "08.25" "09.00"];

jj = 1;
for ii=pump_freq
    data = csvread("test"+i+"/testVIV_"+i+"_"+ii+".csv");
    time = data(:,1);
    disp = data(:,2);
    disp = disp-disp_zero_mean;
    
    [peak peakidx] = findpeaks(disp,'MinPeakDistance',3000,'MinPeakProminence',0.005);
    A(i,jj) = sqrt(2)*rms(disp);
    % log_decrement = log(peak(2)/peak(2+m))/m;
    % zeta(i) = log_decrement/sqrt(4*pi^2+log_decrement^2);
    % f_n(i) = f_d(i)/sqrt(1-zeta(i)^2);
    jj = jj+1;
end
end
% f_d_ave = mean(f_d);
% f_n_ave = mean(f_n);
% f_n_95 = std(f_n)*tinv(0.975,i-1);
% zeta_ave = mean(zeta);
% zeta_95 = std(zeta)*tinv(0.975,i-1);


%% A_star
% U_r = [5 10 15];
D = 0.08;
% U = U_r*f_n_ave*D;
A_star = A/D;
zero_set = zeros(i,1);
A_star = [zero_set A_star];

A_star_ave = mean(A_star,1);
A_star_max = max(A_star,[],1);
A_star_min = min(A_star,[],1);
figure
errorbar([0 str2double(pump_freq)],A_star_ave,abs(A_star_ave-A_star_min),abs(A_star_max-A_star_ave),'.','LineStyle','none','MarkerSize',30,'color','k','LineWidth',3,'MarkerFaceColor','r','MarkerEdgeColor','r')
set(gca, 'FontSize',12)
xlabel('Pump Frequency')
ylabel('A^*')
%% Testing
% drift = csvread("ustar_05hz_test2.csv");
% disp_drift_time = drift(:,1);
% disp_drift_pos = drift(:,2);
% disp_drift_mean = mean(disp_drift_pos);
% figure
% plot(disp_drift_pos)
%% Griffin Plot
% mstar = 7;
% zeta = 0.02;
% C_A = 0.5;
% griffin = (mstar+C_A)*zeta;