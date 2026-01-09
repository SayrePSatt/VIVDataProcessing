%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: error_prop.m
% Version: 1
% Date: 10/23/2025
% Author: Sayre Satterwhite (sayreps@umich.edu)
% Description: Conducts error propagation for appropriate values
% This is adapted from a previous version for use with new data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datafolder = "F:\EFDL\freeDecayTesting\";

%% Specifications
%Directly measured parameters
rho = 998;
d_sph = 0.0889;  %Diameter of Sphere
m_1k = 2.458347-(2/3)*0.0029;
m_6k = 2.458347;
f_s = 1000;     %Sampling Frequency
C_A = 0.5;     %Added mass coefficient
m_d = (4/3)*pi*(d_sph/2)^3*rho+rho*0.005^2*pi*d_sph/4; %Displaced mass

%Indirect measurements
temp_1k = table2array(readtable(datafolder+"1k_09_26_2025/freedecay_1k_air.dat"));
f_n_1k(1,:) = temp_1k(1,:);
f_n_1k_95 = temp_1k(2,1);
% zeta_1k_95 = temp_1k(2,2);
temp_1k = table2array(readtable(datafolder+"1k_09_26_2025/freedecay_1k_water.dat"));
f_w_1k_95 = temp_1k(2,1);
zeta_1k_95 = temp_1k(2,2);
f_w_1k(1,:) = temp_1k(1,:);

temp_6k = table2array(readtable(datafolder+"6k_09_26_2025/freedecay_6k_air.dat"));
f_n_6k(1,:) = temp_6k(1,:);
f_n_6k_95 = temp_6k(2,1);
temp_6k = table2array(readtable(datafolder+"6k_09_26_2025/freedecay_6k_water.dat"));
f_w_6k_95 = temp_6k(2,1);
zeta_6k_95 = temp_6k(2,2);
f_w_6k(1,:) = temp_6k(1,:);

m_a_1k = ((f_n_1k(:,1)./f_w_1k(:,1)).^2-1)*m_1k; %test
m_a_6k = ((f_n_6k(:,1)./f_w_6k(:,1)).^2-1)*m_6k;

St = 0.19;
St_68 = 0.005;
omegana_1k = 2*pi*f_n_1k(:,1);
k_1k = m_1k*omegana_1k.^2; %5.375; %(f_n(1)*2*pi)^2*m
omegana_6k = 2*pi*f_n_6k(:,1);
k_6k = m_6k*omegana_6k.^2;
% k_6k = k_6k - 1.0;
% k_6k(:) = 33.7;

zeta_1k = f_n_1k(:,2);
zeta_6k = f_n_6k(:,2);

c_1k = 4*pi*f_n_1k(:,2).*m_1k.*f_n_1k(:,1);
c_6k = 4*pi*f_n_6k(:,2).*m_6k.*f_n_6k(:,1);
m_star_1k = m_1k/m_d;
m_star_6k = m_6k/m_d;
mass_damp_1k = (m_star_1k+C_A)*f_n_1k(1,2);
mass_damp_6k = (m_star_6k+C_A)*f_n_6k(1,2);
% scruton = 2*m*f_n_1k(2)/(rho*d_sph^2); 

load("pumpFit_freq2velo.mat");

%% Base uncertainties

U_sigma68 = mdl.MSE;
m_sigma68 = 0.5e-3;

f_n_6k_sigma68 = f_n_6k_95/2;
f_n_1k_sigma68 = f_n_1k_95/2;

zeta_6k_sigma68 = zeta_6k_95/2;
zeta_1k_sigma68 = zeta_1k_95/2;

%% Derived uncertainties

k_6k_sigma68 = sqrt((4*pi^2*f_n_6k(1)^2*m_sigma68)^2+(8*pi^2*f_n_6k(1)*m_1k*f_n_6k_sigma68)^2);
k_1k_sigma68 = sqrt((4*pi^2*f_n_1k(1)^2*m_sigma68)^2+(8*pi^2*f_n_1k(1)*m_1k*f_n_1k_sigma68)^2);

c_6k_sigma68 = sqrt((2*sqrt(m_6k*k_6k)*zeta_6k_sigma68)^2 ...
                    +(zeta_6k*sqrt(k_6k/m_6k)*m_sigma68)^2 ...
                    +(zeta_6k*sqrt(m_6k/k_6k)*k_6k_sigma68)^2 );

c_1k_sigma68 = sqrt((2*sqrt(m_1k*k_1k)*zeta_1k_sigma68)^2 ...
                    +(zeta_1k*sqrt(k_1k/m_1k)*m_sigma68)^2 ...
                    +(zeta_1k*sqrt(m_1k/k_1k)*k_1k_sigma68)^2 );