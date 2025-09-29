clear all
close all
clc

m = 1;
fluid = ["water"];
stiffness = 6;

for i=1:length(fluid)
    subfolder = "F:\EFDL\freeDecayTesting\6k_09_26_2025\";
    zero = csvread(subfolder+"freedecay_"+stiffness+"k_"+fluid(i)+"_00.csv");
    disp_zero_time = zero(:,1);
    disp_zero_pos = zero(:,2);
    disp_zero_mean = mean(disp_zero_pos);
    
    for ii=1:1
        data = csvread(subfolder+"freedecay_"+stiffness+"k_"+fluid(i)+"_0"+ii+".csv");
        time = data(:,1);
        disp = data(:,2);
        disp = disp-disp_zero_mean;
        
        [peak peakidx] = findpeaks(disp);
        for jj = 1:length(peakidx)-m
            if time(peakidx(jj)) > 5 && disp(peakidx(jj)) > 0.0005
                f_d(ii,jj) = 1/((time(peakidx(jj+m))-time(peakidx(jj)))/m);
                log_decrement(ii,jj) = log(peak(jj)/peak(jj+m))/m;
            else
                f_d(ii,jj) = NaN;
                log_decrement(ii,jj) = NaN;
            end
        end
        f_d_ave = mean(f_d(ii,:),'omitnan');
        f_d_std =  std(f_d(ii,:),'omitnan');
        zeta(ii,:) = log_decrement(ii,:)./sqrt(4*pi^2+log_decrement(ii,:).^2);
        zeta_ave(ii) = mean(zeta(ii,:),'omitnan');
        zeta_std(ii) =  std(zeta(ii,:),'omitnan');
        f_n(ii,:) = f_d(ii,:)/sqrt(1-zeta_ave(ii)^2);
        f_n_ave(ii) = mean(f_n(ii,:),'omitnan');
        f_n_std(ii) =  std(f_n(ii,:),'omitnan');
        
        peak = peak(1:end-1);
        peakidx = peakidx(1:end-1);
        figure
        hold on
        plot(time,disp,'k')
        plot(time(peakidx),peak,'LineStyle','none','Color','r','Marker','o')
        title("Free Decay Test Time #"+ii)
        figure
        plot(time(peakidx),f_n(ii,:));
        hold on
        yline(f_n_ave(ii),'k');
        yline(f_n_ave(ii)-f_n_std(ii)*2,'r--');
        yline(f_n_ave(ii)+f_n_std(ii)*2,'r--');
        title("Nat. Freq for Test #"+ii);
        figure
        semilogy(time(peakidx),peak,'LineStyle','none','Color','r','Marker','.')
        hold on

    end
    
    f_d_ave = mean(f_d);
    f_n_ave(i) = mean(f_n);
    f_n_95 = std(f_n)*tinv(0.975,ii-1);
    zeta_ave = mean(zeta)
    zeta_95 = std(zeta)*tinv(0.975,ii-1);
    
    % writemtx = [f_n_ave, zeta_ave; f_n_95, zeta_95];
    % writematrix(writemtx,subfolder+"freedecay_"+stiffness+"k_"+type+".dat");
end

%% Reduced velocity calculator

% U_r = [5 10 15];
% D = 0.08;
% U = U_r*f_n_ave*D;

