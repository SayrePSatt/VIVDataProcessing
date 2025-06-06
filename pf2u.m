function [Uinf] = pf2u(pfi) 
%%% Written by Jisheng Zhao, on Mar 11, 2011.
% This is mfile is to calculate the freestream velocity for a given
% pump frequency value, based on interplation of LDV measurement
% data undertaken at a water level of 770 mm on Mar 10, 2010. This
% mfile is mainly written for flow-induced vibration (FIV)
% experimental data processing.
%==============================================

format 'long'
%%% M1: LDV measurement results on 20130311:


%% pf [Hz]      Uinf [m/s]      std(Uinf)       TL (%)
M1=[05.74	0.050308	0.000772	1.533697
   08.26	0.071936	0.001014	1.409894
   10.77	0.093400	0.001286	1.376614
   13.28	0.115380	0.001643	1.424030
   15.78	0.136717	0.001992	1.456916
   18.29	0.158560	0.002267	1.430022
   20.79	0.180274	0.002465	1.367166
   23.29	0.201812	0.002777	1.376147
   25.78	0.222627	0.003167	1.422758
   28.28	0.244162	0.003535	1.447873
   30.77	0.265101	0.003672	1.384997
   33.25	0.286060	0.003752	1.311453
   35.73	0.307336	0.003732	1.214168
   38.21	0.329080	0.004583	1.392737
   40.68	0.350495	0.005456	1.556672
   43.11	0.371510	0.005407	1.455411];
%%% M2: LDV measurement results on 20130313:

M2=[09.69	0.083992	0.001341	1.596965
   11.96	0.103821	0.001690	1.627411
   16.49	0.142871	0.002333	1.632955
   18.76	0.162485	0.002641	1.625524
   21.02	0.181944	0.002967	1.630629
   23.28	0.201356	0.003275	1.626417
   25.54	0.220806	0.003627	1.642622
   27.80	0.240409	0.003894	1.619576
   30.06	0.259891	0.004406	1.695495
   32.30	0.278983	0.004613	1.653636
   34.54	0.298201	0.005032	1.687438
   36.79	0.317575	0.005382	1.694855
   39.03	0.336788	0.005266	1.563520
   41.25	0.356008	0.005793	1.627346
   43.45	0.375212	0.006028	1.606614];


%%% ----------------------------------------
%pft = [M1(:,1); M2(:,1)];
%Uinft = [M1(:,2); M2(:,2)];
pft = M1(:,1);
Uinft = M1(:,2);


%%% Interplation of the free stream velocity 'Uinf'
ppu = interp1(pft,Uinft,'linear','pp');
Uinf = ppval(ppu,pfi); % to be solved
%%% For testing the interplation
%ppf = [5.5:0.01:50];
%puu = interp1(pft,Uinft,ppf,'linear');
%plot(pft,Uinft,'r*',pfi,Uinf,'go',ppf,puu,'b-');


