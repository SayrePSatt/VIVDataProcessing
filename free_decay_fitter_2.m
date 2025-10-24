clear all
close all
clc

m = 1;
fluid = ["air"];
stiffness = 6;

for i=1:length(fluid)
    subfolder = "F:\EFDL\freeDecayTesting\6k_09_26_2025\";
    zero = csvread(subfolder+"freedecay_"+stiffness+"k_"+fluid(i)+"_00.csv");
    disp_zero_time = zero(:,1);
    disp_zero_pos = zero(:,2);
    disp_zero_mean = mean(disp_zero_pos);
    
    for ii=1:4
        clear f_n f_d peakidx peak log_decrement zeta
        data = csvread(subfolder+"freedecay_"+stiffness+"k_"+fluid(i)+"_0"+ii+".csv");
        time = data(:,1);
        disp = data(:,2);
        disp = disp-disp_zero_mean;
        
        [peak peakidx] = findpeaks(disp,'MinPeakHeight',0);
        for jj = 1:length(peakidx)-m
            if time(peakidx(jj)) > 50 && time(peakidx(jj)) < 200 && disp(peakidx(jj)) > 0.0005
                f_d(jj) = 1/((time(peakidx(jj+m))-time(peakidx(jj)))/m);
                log_decrement(jj) = log(peak(jj)/peak(jj+m))/m;
            else
                f_d(jj) = NaN;
                log_decrement(jj) = NaN;
            end
        end
        f_d_ave(i,ii) = mean(f_d,'omitnan');
        f_d_std(i,ii) =  std(f_d,'omitnan');
        zeta = log_decrement./sqrt(4*pi^2+log_decrement.^2);
        zeta_ave(i,ii) = mean(zeta,'omitnan');
        zeta_std(i,ii) =  std(zeta,'omitnan');
        f_n = f_d./sqrt(1-zeta.^2);
        f_n_ave(i,ii) = mean(f_n,'omitnan');
        f_n_std(i,ii) =  std(f_n,'omitnan');
        
        peak = peak(1:end-1);
        peakidx = peakidx(1:end-1);
        
        figure
        subplot(2,2,1)
        hold on
        plot(time,disp,'k')
        plot(time(peakidx),peak,'LineStyle','none','Color','r','Marker','o')
        title("Free Decay Test Time #"+ii)

        subplot(2,2,2)
        plot(time(peakidx),f_n);
        hold on
        yline(f_n_ave(i,ii),'k');
        yline(f_n_ave(i,ii)-f_n_std(i,ii)*2,'r--');
        yline(f_n_ave(i,ii)+f_n_std(i,ii)*2,'r--');
        title("Nat. Freq for Test #"+ii);
        
        subplot(2,2,4)
        plot(time(peakidx),zeta);
        hold on
        yline(zeta_ave(i,ii),'k');
        yline(zeta_ave(i,ii)-zeta_std(i,ii)*2,'r--');
        yline(zeta_ave(i,ii)+zeta_std(i,ii)*2,'r--');
        title("Damp. Ratio for Test #"+ii);

        subplot(2,2,3)
        semilogy(time(peakidx),peak,'LineStyle','none','Color','r','Marker','.')
        hold on

    end
    
    f_d_ave_all(i) = mean(f_d_ave(i,:));
    f_n_ave_all(i) = mean(f_n_ave(i,:));
    f_n_95(i) = mean(f_n_std(i,:))*2;
    zeta_ave_all(i) = mean(zeta_ave(i,:));
    zeta_95(i) = mean(zeta_std(i,:))*2;
    
    writemtx = [f_n_ave_all(i), zeta_ave_all(i); f_n_95(i), zeta_95(i)];
    writematrix(writemtx,subfolder+"freedecay_"+stiffness+"k_"+fluid(i)+".dat");
end

%% Reduced velocity calculator

% U_r = [5 10 15];
% D = 0.08;
% U = U_r*f_n_ave*D;

