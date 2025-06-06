function [phase_TT,phase_VT]= retrievephase(yN,CyT,CvT)

  Hyy = hilbert(yN);
  %HCyM = hilbert(CyM);
  HCyT = hilbert(CyT);
  %HCvM = hilbert(CvM);
  HCvT = hilbert(CvT);
  %% total phase and vortex phase
  ang_yy = angle(Hyy)*360/(2*pi);
  %ang_CyM = angle(HCyM)*360/(2*pi);
  ang_CyT = angle(HCyT)*360/(2*pi);
  %ang_CvM = angle(HCvM)*360/(2*pi);
  ang_CvT = angle(HCvT)*360/(2*pi);
  %%
  %phase_TM = ang_CyM - ang_yy';
  phase_TT = ang_CyT - ang_yy;
  %phase_VM = ang_CvM - ang_yy';
  phase_VT = ang_CvT - ang_yy;
%   for N = 1:length(phase_TM)
%     if 270<=phase_TM(N)
%       phase_TM(N)=phase_TM(N)-360;
%     elseif phase_TM(N)<=-90
%       phase_TM(N)=phase_TM(N) + 360;
%     end
%   end
   for N = 1:length(phase_TT)
    if 270<=phase_TT(N)
      phase_TT(N)=phase_TT(N)-360;
    elseif phase_TT(N)<=-90
      phase_TT(N)=phase_TT(N) + 360;
    end
  end
  
%   for N = 1:length(phase_VM)
%     if 270<=phase_VM(N)
%       phase_VM(N)=phase_VM(N)-360;
%     elseif phase_VM(N)<=-90
%       phase_VM(N)=phase_VM(N) + 360;
%     end
%   end
  
  
  for N = 1:length(phase_VT)
    if 270<=phase_VT(N)
      phase_VT(N)=phase_VT(N)-360;
    elseif phase_VT(N)<=-90
      phase_VT(N)=phase_VT(N) + 360;
    end
  end
  
  
  
