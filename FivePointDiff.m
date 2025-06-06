function dx=FivePointDiff(x,fs)
% 
% fs is the sampling frequency
h = 1/fs;
for i=1:length(x)-5
dx(i)=(-25*x(i)+48*x(i+1)-36*x(i+2)+16*x(i+3)-3*x(i+4))/(12*h);
end