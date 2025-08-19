clear all
close all
clc

load("pumpFit_velo2freq.mat")
stiffness = 6;

temp = table2array(readtable("E:\EFDL\tandemSphere\freeDecay/"+stiffness+"k_08_18_2025/freedecay_"+stiffness+"k_water.dat"));
f_w = temp(1,1);

d = 0.0889;
if stiffness==1
    U_star = [6.5:1:23.5]';
elseif stiffness==6
    U_star = [2.5:0.5:6.5]';
end

U_star = flip(U_star);
U = U_star*(f_w*d);

freq = predict(mdl,U);

transient_duration = 200;
transient_duration = transient_duration*ones(size(freq));

txtfile = [0 10;freq transient_duration];
txtfile = round(txtfile,2);

writematrix(txtfile,"E:/EFDL/test_specs/"+stiffness+"k_tests_0.5_reverse_08_18_2025.txt","Delimiter",',')
% save("E:/EFDL/test_specs/6k_tests_0.5_08_18_2025.txt","txtfile","-ascii","-tabs")