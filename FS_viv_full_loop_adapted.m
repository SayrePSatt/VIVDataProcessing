%% CODE FOR PROCESSING VIV OF A SPHERE

clear all;
close all;

addpath('E:/EFDL/Old_Data');
mypath = 'E:/EFDL/Old_Data';

%% immersed lengths to process 
%IL = [00.000 00.125 00.250 00.375 00.500 01.000 01.250 01.500 01.750 02.000 ];

IL = 01.00;
for k = 1:length(IL)
%% Measuring reference value (taring the linear encoder) 

 
  B = importdata([mypath sprintf('/RawData/001_pf00.00hz.dat', IL(k))], '\t',23); dat0 = B.data;
  lvdt0 = mean(dat0(:,5));
    
    
%% reading fna,zeta and fnw from the file
path1= strcat([mypath sprintf('/freedecay_air/freedecay.dat', IL(k))]);
fileID = fopen(path1, 'r');
fn1 = textscan(fileID, '%f%f');
fna = fn1{1}
zetaair= fn1{2}
fclose(fileID);

path2= strcat([mypath sprintf('/freedecay_water/freedecay.dat', IL(k))]);
fileID = fopen(path2, 'r');
fn2 = textscan(fileID, '%f%f');
fnw = fn2{1}; fnw = fnw(1);
zetawater= fn2{2};zetawater = zetawater(1);
fclose(fileID);
%%  reading total mass and mass bottom from the dat file
% path3= strcat([mypath '\mass.dat']);
% fid = fopen(path3, 'r');
% m =textscan(fid, '%f%f');
% mass = m{1};
% mbt = m{2};
mass =   1392.8/1000;  %%% NOTE: MASS TO BE CALCULATED BY THE END OF THE EXPERIMENT
mbt = 345.7/1000;
%% DEFINING PARAMETERS OF THE PROBLEM


DD = 0.08; % diameter of the sphere [m]
rho = 998; % density of water [kg/m3]
massd = 4/3*pi*(DD/2)^3*rho; % displaced mass of the sphere [kg]
omegana=2*pi*fna;
m_star= mass/massd;
c_damping = 2*mass*zetaair*omegana; % structural damping or damping coefficient in air (c)
kspr=mass*omegana^2;% spring constant measured in air
mass_A = ((fna/fnw)^2 - 1)*mass;% theoretical added mass
mtp=mass-mbt; % top mass
kl = 1; % Linear encoder
Fs = 100; %sampling frequency
CA = 0.5; % added mass coefficient is 1 for cylinder
% ---------------------------------------------------------------------------------------

datfiles = dir([mypath sprintf('/RawData/*pf*.dat', IL(k))]);
j = 0;
for i = 2:length(datfiles)-1
     name = datfiles(i).name;
    j = j+1;
    if(findstr(name, 'pf00.00'))  % skipping the reference files
        checkname = 1;
    else
        checkname =0;
    end
    if(checkname==1)
        continue
    end
    
    % Extracting linear encoder data and subtracting the reference value
    A = importdata([mypath sprintf('/RawData/', IL(k)) name], '\t',23); dat = A.data;
    t = dat(:,1); % time
    lvdt_raw = dat(:,5); % lvdt signal
    lvdt_raw(isnan(lvdt_raw))=0;
    
    mm = (lvdt_raw-lvdt0)*kl; % referencing the LVDT signal
    
    % extracting pump frequency from the name and calculating ustar
    pf = str2num(name(end-10:end-6));
    u(j) = pf2u(pf);
    ustar(j) = u(j)./(fnw*DD);
    
    Pstat = 0.5*998*u(j)^2*pi*(DD/2)^2; % static pressure
    
    %% BUTTERWORTH FILTERING
    
    Fs = 100; % samping freq.
    Fc = 2; % [Hz], cutoff freq. of Butterworth filter
    Wn_cut = 2*Fc/Fs; %The cutoff frequency Wn must be 0.0 < Wn < 1.0 with 1.0 corresponding to half the sample rate.
    N = 4; % Nth order lowpass digital
    
    [B,A] = butter(N,2*Fc/Fs,'low'); % 4-order Butterworth
    umm=mm; % unfiltered mm
    mm = filtfilt(B,A,mm);% filtered mm
    
    %mm1 = mm;
    
    %compare filtered and unfiltered to check if you are filtering right
%   plot(t, mm, 'b', t, umm, 'k'); xlim([200 220]);
%   pause
    
    
    
    ymean(j) = mean(mm);
    mm = mm -ymean(j); % subtract the mean of the displacement

    y1 = mm/(DD*1000); % non-dimensional displacement
    y1 = y1(1:end-10);
    
    %% FINDING Y_RMS AS GIVEN IN GOVARDHAN 2005
    
    yrms= rms(mm); % take the rms of the fluctuations
    ymax = max(abs(mm)); 
    
    P(j)= sqrt(2)*yrms./ymax;
    
    Arms(j)= sqrt(2)* yrms./(DD*1000); % definition of A_star used by Govardhan and Williamsons
    Amean(j) = ymean(j)/(DD*1000);
%----------------------------------------
    % mean of top 10% (Some studies use this definition)
     
%         [Ax Bx]=peakdet((mm1-mean(mm1)),std(mm1));
%         
%         [ap kp]=sort(Ax(:,2),'descend'); % sort Ax(+) from highest
%         [an kn]=sort(-Bx(:,2),'descend'); % sort -Bx(-) from highest
%         kp=Ax(kp,1);
%         kn=Bx(kn,1);        	
%         wwp = ceil(length(kp)*0.10); % Top ten percent of kp
%         wwn = ceil(length(kn)*0.10);
%         a10n=mean([mm1(kn(1:wwn))]); % Top ten percent indices in wwn to top ten percent of values in yf, then mean. Mean of the top ten percent
%         a10p=mean([mm1(kp(1:wwp))]);
%         a10=mean([mm1(kp(wwp)) - mm1(kn(wwn))]); % Top ten percent (+) - top ten percent (-)
%         Amax10(j) = a10p./(DD*1000);
%         clearvars a kk Ax Bx ap kp an kn kp wwp wwn a10n a10p a10   
%            
%----------------------------------------      
    
    %% FREQUENCY ANALYSIS
    
    [yypwr yyf] = psdd3(Fs,mm,2048,1024,1);
    yyf1 = yyf;
    yyf = yyf ./fnw;% normalized frequency
    
    temp_meanpwr = mean(yypwr,2)';
    [yypwr_mx yyf_mx] = max(temp_meanpwr);
    
    p_yy = polyfit(yyf(yyf_mx-1:yyf_mx+1),yypwr(yyf_mx-1:yyf_mx+1),2);
    
    fpos = -p_yy(2)/(2*p_yy(1)); %1st derivative=0, find freq value
    
    p_yy1 = polyfit(yyf1(yyf_mx-1:yyf_mx+1),yypwr(yyf_mx-1:yyf_mx+1),2);
    
    fpos1 = -p_yy1(2)/(2*p_yy1(1));
    
    fstar(j) = fpos;
    
    
    T(j) = 1/fpos;% time period of oscillation
    
    t_cycl = t/T(j);
    t_cycl =t_cycl(1:end-10); % for getting the time trace of y/D
     
    % plotting the frequency analysis: 
    %semilogy(yyf, yypwr); xlim([0 4]);  
    %% TIME TRACE OF THE SIGNAL 
    plot(t_cycl, y1); xlim([0 20]); title(sprintf('$U^*$ = %2.1f', ustar(j)));% time trace of the displacement
     %pause
     %% CALCULATING EFFECTIVE ADDED MASS FOR THE SPHERE:
    
    C_EA(j) = m_star*((1-fstar(j).^2)/fstar(j).^2) + (CA/fstar(j).^2); 
    
    %% MASS DAMPING OF THE SYSTEM:
    CA_2 = ((fna/fnw)^2 - 1)*m_star;
    Massdamping(j) = (m_star+CA_2)*zetaair;
    
    %% CALCULATING ACCELERATION FROM THE DISPLACEMENT
    
    y = mm/1000; % displacement in metres
    tt = t(1:end-10);
    yv = FivePointDiff(y,Fs);
    yv = filtfilt(B,A,yv);%%% smoothing
    
    ya = FivePointDiff(yv,Fs);
    ya = filtfilt(B,A,ya);
    
    
    yv(end-4:end)=[];
    y(end-9:end)=[];
    
    yN = mm/(DD*1000); % Normalised by the sphere's diameter
    yN2 = yN(1:end-10); % removing last ten points
    

    
    %% THEORETICAL LIFT FORCES CALCULATIONS
    kcor = 1;
    
    FyT = mass*ya + c_damping*yv + kspr*kcor*y'; % theoretical force from theory (it is still derived from measured y, mass, c_damping and kspring
    FyT = filtfilt(B,A,FyT); % filtering the theoretical force\
    CyT = FyT/Pstat;
   
    
    %% COMPARE MEAN FOR ROTATING SPHERE AND RMS FOR NON ROTATING SPHERE
    
  
    % non-rotating sphere
    
    FyT_rms(j) = rms(FyT);
    CyT_rms(j) = FyT_rms(j)/Pstat;
    FyT_mean(j) = mean(FyT);
    CyT_mean(j) = FyT_mean(j)/Pstat;
      %
    %% CALCULATING MEASURED AND THEORETICAL POTENTIAL FORCES
    FpotM= -mass_A *ya; % measured potential forces
    FpotT = -CA * massd * ya; % theoretical potential forces
  
    
    %%  CALCULATING THEORETICAL AND MEASURED VORTEX FORCES
    FvT = FyT-FpotT; %theoretical vortex force when we subtract Fpot from thoretical total force
     % measured vortex force
    CvT = FvT /Pstat;
    
    FvT_rms(j) = rms(FyT);
  
    CvT_rms(j) = FvT_rms(j)/Pstat;
   
    
    %% PLOTTING THE TIME TRACE (2 CYCLES) FOR ALL THE COEFFICIENTS. USE PAUSE WITH THIS TO GET THE RIGHT USTAR
% 
%         figure(5)
%         subplot(4,1,1)
%         plot(t_cycl, y1); ylabel('$y/D$');title(sprintf('$U^* = %3.1f$',ustar(j))); xlim([20 22]); 
%         subplot(4,1,2)
%         plot(t_cycl, CyM); ylabel('$C_{total}$'); xlim([20 22]); 
%         subplot(4,1,3)
%         plot(t_cycl,CpotM); ylabel('$C_{potential}$'); xlim([20 22]);
%         subplot(4,1,4)
%         plot(t_cycl, CvM); ylabel('$C_{vortex}$');xlim([20 22]);xlabel('$t/T$'); 
    
    
    
    
    
    %    plot(t_cycl, y1); xlim([0 20]); title(sprintf('$U^* = %3.1f$',ustar(j)));ylabel('$y/D$');xlabel('$t/T$')



    
    
    %% DEFINING (U*/f*)S
    ustar_s(j) = (ustar(j)./fstar(j))*0.175; % S = 0.2 for the sphere at this Reynomassd = 4/3*pi*(DD/2)^3*rholds number
    

    %% CALCULATING PHASE ANGLES
    
   [phase_TT,phase_VT]= retrievephase1(yN2,CyT,CvT);	
   
    mphi_TT(j) = mean(phase_TT);       
          
    		        
    mphi_VT(j) = mean(phase_VT);
    
%   pause
    
end



%% PLOTTING


   
plotStyle = {'bo--', 'rs--', 'gd--','kv--', 'y^--', 'c>--', 'm<--', 'p--', '+--', 'd--'};
plotStyle1 = {'bo', 'rs', 'gd','kv', 'y^', 'c>', 'm<', 'p'};
legendInfo{k} = sprintf('$IL$ : %4.3f', IL(k)) ; % or whatever is appropriate


figure (3)
h3 = plot(ustar, Arms, plotStyle{k}, 'Linewidth', 2);% set(h3, 'MarkerFaceColor', get(h3, 'Color'));
ylabel('$A_{rms}$');xlabel('$U^*$');set(get(gca, 'YLabel'), 'Rotation', 0);
hold on;

figure(4)% fstar
h4 = plot(ustar, fstar, plotStyle{k}, 'Linewidth', 2); %set(h4, 'MarkerFaceColor', get(h4, 'Color'));
ylabel('$f^*$');xlabel('$U^*$');set(get(gca, 'YLabel'), 'Rotation', 0);
hold on;
% 
figure(5)% fstar
h5 = plot(ustar, Amean, plotStyle{k}, 'Linewidth', 2);%set(h5, 'MarkerFaceColor', get(h5, 'Color'));
ylabel('$\overline{y}/D$');xlabel('$U^*$');set(get(gca, 'YLabel'), 'Rotation', 0);
hold on;

figure(6)
h6 = plot(ustar, P, plotStyle{k}, 'Linewidth', 2); %set(h6, 'MarkerFaceColor', get(h6, 'Color'));
xlabel('$U^*$'); ylabel('Periodicity');set(get(gca, 'YLabel'), 'Rotation', 0);
hold on;

figure(7)
h7 = plot(ustar, CyT_mean, plotStyle{k}, 'Linewidth', 2); 
ylabel('$\overline{C}_y$');xlabel('$U^*$');set(get(gca, 'YLabel'), 'Rotation', 0);
 hold on;

figure(8)
h8 = plot(ustar, CyT_rms, plotStyle{k}, 'Linewidth', 2); 
ylabel('${C_y}_{rms}$');xlabel('$U^*$');set(get(gca, 'YLabel'), 'Rotation', 0);
title('Theoretical'); hold on;


figure(9)
h9 = plot(ustar, mphi_TT, plotStyle{k}, 'Linewidth', 2);
ylabel('$\phi_{total}$');xlabel('$U^*$');set(get(gca, 'YLabel'), 'Rotation', 0);title('Theoretical');
hold on;

figure(10)
h10 = plot(ustar, mphi_VT, plotStyle{k}, 'Linewidth', 2);
ylabel('$\phi_{vortex}$');xlabel('$U^*$');set(get(gca, 'YLabel'), 'Rotation', 0);title('Theoretical');
hold on;


figure(11)
h11 = plot(ustar, C_EA, plotStyle{k}, 'Linewidth', 2);ylabel('$C_{EA}$'); xlabel('$U^*$');
set(get(gca, 'YLabel'), 'Rotation', 0);hold on;


%CLEAR ALL THE VARIABLES BEFORE IT MOVES TO THE NEXT ALPHA LOOP OTHERWISE
%IT WILL ONLY OVERWRITE AND ASSUME THE PREVIOUS VALUES FOR THE OTHER POINTS
%VERY VERY IMPORTANT TO CLEAR VARIABLES (clear all the variables changing
%values with the alpha loop

clearvars Arms fstar Amean P CyT_rms   ustar mphi_TT  mphi_VT C_EA  CyT_mean  mphi_TT  mphi_VT 




end   


figure(3); legend(legendInfo);
saveas(figure(3),[mypath '/figures/comp_Arms.fig']);
print([mypath '/figures/comp_Arms'], '-deps', '-r300');

figure(4); legend(legendInfo);
saveas(figure(4),[mypath '/figures/comp_f.fig']);
print([mypath '/figures/comp_f'], '-deps', '-r300');

figure(5); legend(legendInfo);
saveas(figure(5),[mypath '/figures/comp_Amean.fig']);
print([mypath '/figures/comp_Amean'], '-deps', '-r300');

figure(6); legend(legendInfo)
saveas(figure(6),[mypath '/figures/comp_P.fig']);
print([mypath '/figures/comp_P'], '-deps', '-r300');

figure(7); legend(legendInfo);
saveas(figure(7),[mypath '/figures/comp_CyTmean.fig']);
print([mypath '/figures/comp_CyTmean'], '-deps', '-r300');

figure(8); legend(legendInfo);
saveas(figure(8),[mypath '/figures/comp_CyTrms.fig']);
print([mypath '/figures/comp_CyTrms'], '-deps', '-r300');


figure(9); legend(legendInfo);
saveas(figure(9),[mypath '/figures/comp_phi_TT.fig']);
print([mypath '/figures/comp_phi_TT'], '-deps', '-r300');

figure(10); legend(legendInfo);
saveas(figure(10),[mypath '/figures/comp_phi_VT.fig']);
print([mypath '/figures/comp_phi_VT'], '-deps', '-r300');

figure(11); legend(legendInfo);
saveas(figure(11),[mypath '/figures/comp_CEA.fig']);
print([mypath '/figures/comp_CEA'], '-deps', '-r300');

