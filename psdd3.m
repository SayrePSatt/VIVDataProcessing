function [mx,f] = psdd3(fs,x,N,Ns,op)
% [mx,f] = PSDD3_JV(fs,x,N,Ns,op)
% Windowed FFT.
% fs - sampling frequency of original signal
% x  - time series
% N  - window length
% Ns - window shift
% op - [1] is for Hanning method (default), [2] is Hamming
%
% PSDD3 returns the frequency information mx from each block of data taken from
% x. The code first splits the original time series up into N-long pieces
% with a shift of Ns. Then, windowing is performed on each block, and an
% FFT is calculated and returned in fs.
%
% History:
% originally written by David
% updated April 2013 by James using a new method to break up the data.
%
% See also PSDD

for ii = 1:ceil((length(x)-N+1)/Ns);                % loop for each window
    xx(:,ii) = x(((ii-1)*Ns+1):((ii-1)*Ns+N));      % current window's x data
end

% Windowing:
if op == 1
    YY=(hanning(N)*ones(1,size(xx,2))).*(xx-ones(N,1)*(mean(xx)));
else
    YY=(hamming(N)*ones(1,size(xx,2))).*(xx-ones(N,1)*(mean(xx)));
end

nfft = 2^(nextpow2(N)+2);

% FFTing each column of xx
fftx = fft(YY,nfft);
NumUniquePts = ceil((nfft+1)/2);
fftx = fftx(1:NumUniquePts,:);
mx = abs(fftx).^2/N/fs;
f = (0:NumUniquePts-1)*fs/nfft;
f=f'*ones(1,size(mx,2));
%Tx=ones(size(mx,1),1)*TI;
end