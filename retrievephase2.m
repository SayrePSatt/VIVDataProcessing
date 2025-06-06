function [phase_TT,phase_VT]= retrievephase(fs,freq,yN,CyT,CvT)

[crossCorr_T, lags_T] = xcorr(CyT, yN);
[crossCorr_V, lags_V] = xcorr(CvT, yN);

% Find the index of the maximum cross-correlation
[~, maxIndex_T] = max(crossCorr_T);
[~, maxIndex_V] = max(crossCorr_V);

% Corresponding lag (in samples) between the signals
timeLag_T = lags_T(maxIndex_T);
timeLag_V = lags_T(maxIndex_V);


% Compute the phase lag in radians
phaseLagRad_T = 2 * pi * freq * timeLag_T / fs;
phaseLagRad_V = 2 * pi * freq * timeLag_V / fs;

% Convert phase lag to degrees
phase_TT = rad2deg(phaseLagRad_T);
phase_VT = rad2deg(phaseLagRad_V);

phase_TT = mod(phase_TT, 360);  % Get result in the range [0, 360)
if (phase_TT > 180)
    phase_TT = 360 - phase_TT;
end

phase_VT = mod(phase_VT, 360);  % Get result in the range [0, 360)
if (phase_VT > 180)
    phase_VT = 360 - phase_VT;
end
% 
% for N = 1:length(phase_VT)
%     if 270<=phase_VT(N)
%       phase_VT(N)=phase_VT(N)-360;
%     elseif phase_VT(N)<=-90
%       phase_VT(N)=phase_VT(N) + 360;
%     end
% end
% 
% for N = 1:length(phase_TT)
%     if 270<=phase_TT(N)
%       phase_TT(N)=phase_TT(N)-360;
%     elseif phase_TT(N)<=-90
%       phase_TT(N)=phase_TT(N) + 360;
%     end
% end