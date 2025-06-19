clear all
close all
clc

D = 0.0889;
f_w = 0.23931377;

u_star = [0 linspace(6.5,23.5,35)]

for i=1:length(u_star)
    pumpFreq(i) = pumpSpeedCalculator(u_star(i)*D*f_w);
end

pumpFreq'