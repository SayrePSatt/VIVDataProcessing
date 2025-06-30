clear all
close all
clc

D = 0.0889;
f_w = 0.56804767708258;

u_star = [0 2.5:2:6.5]

for i=1:length(u_star)
    pumpFreq(i) = pumpSpeedCalculator(u_star(i)*D*f_w);
end

pumpFreq'