clear all
close all
clc

m = 5;
type = "air";
subfolder = "free_decay_90mm_1k_alt\";
zero = csvread(subfolder+"freedecay_"+type+"_00.csv");
disp_zero_time = zero(:,1);
disp_zero_pos = zero(:,2);
disp_zero_mean = mean(disp_zero_pos);

for ii=1:8
    data = csvread(subfolder+"freedecay_"+type+"_0"+ii+".csv");
    time = data(:,1);
    disp = data(:,2);
    disp = disp-disp_zero_mean;
    
    [peak peakidx] = findpeaks(disp);
    f_d(ii) = 1/((time(peakidx(2+m))-time(peakidx(2)))/m);
    log_decrement = log(peak(2)/peak(2+m))/m;
    zeta(ii) = log_decrement/sqrt(4*pi^2+log_decrement^2);
    f_n(ii) = f_d(ii)/sqrt(1-zeta(ii)^2);
end

f_d_ave = mean(f_d);
f_n_ave = mean(f_n);
f_n_95 = std(f_n)*tinv(0.975,ii-1);
zeta_ave = mean(zeta);
zeta_95 = std(zeta)*tinv(0.975,ii-1);

writemtx = [f_n_ave, zeta_ave; f_n_95, zeta_95];
writematrix(writemtx,subfolder+"freedecay_"+type+".dat");

findpeaks(disp)

scrouton = 2*2*zeta_ave/(1000*0.08^2);

%% Reduced velocity calculator
U_r = [5 10 15];
D = 0.08;
U = U_r*f_n_ave*D;