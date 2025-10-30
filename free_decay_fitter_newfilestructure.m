clear all
close all
clc

m = 1;
fluid = ["air"];
stiffness = ["6k"];
folder = uigetdir;
allfiles = dir(folder);
[~, date, ~] = fileparts(folder);
date = string(date);
%%
idx = 1;
for i=3:length(allfiles);
    filename = "\"+string(allfiles(i).name);
    clear f_d f_n log_decrement
    if contains(filename,fluid) && endsWith(filename,'.dat') && contains(filename,stiffness) && ~contains(filename,"_0.0_")
        data = readtable(folder+filename);
        data = table2array(data);
        time = data(:,1);
        disp = data(:,2);
        [peak, peakidx] = findpeaks(disp,'MinPeakHeight',0);
        for jj = 1:length(peakidx)-m
            if time(peakidx(jj)) > 30 && time(peakidx(jj)) < 200 && disp(peakidx(jj)) > 0.0005
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
    end
end
f_d_ave_all = std_multimean(num_samples,f_d_ave,f_d_std);
f_n_ave_all = std_multimean(num_samples,f_n_ave,f_n_std);
zeta_ave_all = std_multimean(num_samples,zeta_ave,zeta_std);

T = table(f_n_ave_all,zeta_ave_all,f_d_ave_all,'VariableNames',{'f_n', 'zeta', 'f_d'});

% writemtx = [f_n_ave_all(idx), zeta_ave_all(idx); f_n_95(idx), zeta_95(idx)];
writetable(T,folder+"\"+date+"_freedecay_"+stiffness+"_"+fluid(1)+".dat");

%% Reduced velocity calculator

% U_r = [5 10 15];
% D = 0.08;
% U = U_r*f_n_ave*D;

