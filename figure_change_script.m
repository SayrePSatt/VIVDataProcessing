function figure_change_script(change_fig,font_size,mark_size)
% Defaults for this blog post
%addpath('/nas-6000/wclab/anchal/octave-lib/');
%cd /nas-6000/wclab/anchal/
%addpath('/nas-6000/wclab/anchal/matlab-lib');

for ii = 1:length(change_fig)
    figure(change_fig(ii))

% Defaults for this blog post
width = 5;     % Width in inches
height = 4.5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = font_size;      % Fontsize
lw = 2.0;      % LineWidth
msz = 8;       % MarkerSize

h = findobj(gca, 'Type', 'Line');
set(h, 'MarkerSize', mark_size)

set(gca, 'FontSize', fsz);
set(0,'defaulttextinterpreter','latex')

% % The properties we've been using in the figures
% set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
% set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
% set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
% set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
% 
% % Set the default Size for display
 defpos = get(0,'defaultFigurePosition');
 set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*120, height*100]);
end
% 
% % Set the defaults for saving/printing to a file
% set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
% set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
% defsize = get(gcf, 'PaperSize');
% left = (defsize(1)- width)/2;
% bottom = (defsize(2)- height)/2;
% defsize = [left, bottom, width, height];
% set(0, 'defaultFigurePaperPosition', defsize);
end