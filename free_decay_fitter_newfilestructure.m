clear all
close all
clc

m = 1;
stiffness = ["1k"];
folder = uigetdir;
allfiles = dir(folder);
[~, date, ~] = fileparts(folder);
date = string(date);
%%
fluid = ["water"];
idx = 1;
for i=3:length(allfiles)
    filename = "\"+string(allfiles(i).name)
    clear f_d f_n log_decrement data
    searchbool = contains(filename,fluid) && endsWith(filename,'.dat') && contains(filename,stiffness) && ~contains(filename,"results");
    if  searchbool && ~contains(filename,"_00.0_")
        data = readtable(folder+filename);
        data = table2array(data);
        time = data(:,1);
        disp = data(:,2);
        [peak, peakidx] = findpeaks(disp,'MinPeakHeight',0,'MinPeakDistance',300);
        for jj = 1:length(peakidx)-m
            if time(peakidx(jj)) > 50 && time(peakidx(jj)) < 200 && disp(peakidx(jj)) > 0.0005
                f_d(jj) = 1/((time(peakidx(jj+m))-time(peakidx(jj)))/m);
                log_decrement(jj) = log(peak(jj)/peak(jj+m))/m;
            else
                f_d(jj) = NaN;
                log_decrement(jj) = NaN;
            end
        end
        num_samples(idx) = length(f_d(~isnan(f_d)));
        f_d_ave(idx) = mean(f_d,'omitnan');
        f_d_std(idx) =  std(f_d,'omitnan');
        zeta = log_decrement./sqrt(4*pi^2+log_decrement.^2);
        zeta_ave(idx) = mean(zeta,'omitnan');
        zeta_std(idx) =  std(zeta,'omitnan');
        f_n = f_d./sqrt(1-zeta.^2);
        f_n_ave(idx) = mean(f_n,'omitnan');
        f_n_std(idx) =  std(f_n,'omitnan');

        peak = peak(1:end-1);
        peakidx = peakidx(1:end-1);

        figure
        subplot(2,2,1)
        hold on
        plot(time,disp,'k')
        plot(time(peakidx),peak,'LineStyle','none','Color','r','Marker','o')
        title("Free Decay Test Time #"+i)

        subplot(2,2,2)
        plot(time(peakidx),f_n);
        hold on
        yline(f_n_ave(idx),'k');
        yline(f_n_ave(idx)-f_n_std(idx)*2,'r--');
        yline(f_n_ave(idx)+f_n_std(idx)*2,'r--');
        title("Nat. Freq for Test #"+idx);

        subplot(2,2,4)
        plot(time(peakidx),zeta);
        hold on
        yline(zeta_ave(idx),'k');
        yline(zeta_ave(idx)-zeta_std(idx)*2,'r--');
        yline(zeta_ave(idx)+zeta_std(idx)*2,'r--');
        title("Damp. Ratio for Test #"+idx);

        subplot(2,2,3)
        semilogy(time(peakidx),peak,'LineStyle','none','Color','r','Marker','.')
        hold on
        idx = idx + 1;
    elseif searchbool && contains(filename,"_00.0_")
        data = readtable(folder+filename,'NumHeaderLines',9,'Range','A:B');
        data = table2array(data);
        data = data(1:2,:);
        mass = data(:,1);
        diameter = data(:,2);
    end
end
f_d_ave_all = std_multimean(num_samples,f_d_ave,f_d_std);
f_n_ave_all = std_multimean(num_samples,f_n_ave,f_n_std);
zeta_ave_all = std_multimean(num_samples,zeta_ave,zeta_std);

T = table(f_n_ave_all,zeta_ave_all,mass,diameter,f_d_ave_all,'VariableNames',{'f_n (Hz)', 'zeta', 'mass (kg)', 'diameter (m)', 'f_d (Hz)'});

% writemtx = [f_n_ave_all(idx), zeta_ave_all(idx); f_n_95(idx), zeta_95(idx)];
writetable(T,folder+"\"+date+"_freedecay_results_"+stiffness+"_"+fluid(1)+".dat");

%% Reduced velocity calculator

% U_r = [5 10 15];
% D = 0.08;
% U = U_r*f_n_ave*D;

