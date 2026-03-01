clear all
close all
clc

mkr_sz = 8;
line_wdth = 1.2;
fig_sz = [100 100 600 350];
tick_sz = [0.03 0.012];
font_sz = 20;
lgd = 0;
xlbl = 0;
ylbl = 0;
titl = 0;
xscl = 1;
yscl = 1;

sub_fig = 1;
rows = 

load("pumpFit_freq2velo.mat");

topfolder = 'C:\Users\sayre\Documents\Github\VIVProcessing\figures\matlab_figs\';
if sub_fig == 1
    [processing_files, processing_file_dir] = uigetfile(topfolder,'*.fig','MultiSelect','on');
else
    [processing_files, processing_file_dir] = uigetfile(topfolder,'*.fig','MultiSelect','on');
end


for ii = 1:length(processing_files)
    close all
    current_fig = openfig(fullfile(processing_file_dir,processing_files{ii}));
    figure(current_fig)
    hold on
    ax = gca;
    ax.ActivePositionProperty = 'position';
    pos = ax.Position; % Save current axes position
    
    if lgd == 0
        lgd = findobj(gcf, 'Type', 'Legend');  % Find legend in the figure
        lgd.Visible = 'off';
    end

    ax.Position = pos
    if xlbl == 0
        xlabel(' ')
        gca.Position = pos;
    end

    if ylbl == 0
        ylabel(' ')
        gca.Position = pos;
    end

    % if titl == 0
    %     titl(' ')
    %     gca.Position = pos;
    % end
    
    if xscl == 0
        xticklabels({' '})
        gca.Position = pos;
    end

    if yscl == 0
        yticklabels({' '})
        gca.Position = pos;
    end

    h = findobj(current_fig, 'Type', 'Line');   % This finds all line objects (including markers)
    
    % Change the marker size for all found lines
    for k = 1:numel(h)
        h(k).MarkerSize = mkr_sz;
        h(k).LineWidth = line_wdth;
    end

    ax = findobj(current_fig,'Type','Axes');
    set(ax,'FontSize',font_sz)
    set(ax,'LineWidth',line_wdth)

     current_fig.Position = fig_sz;

    exportgraphics(current_fig,['figures\' processing_files{ii}(1:end-4) '.pdf'],'Resolution',300)
    exportgraphics(current_fig,['figures\' processing_files{ii}(1:end-4) '.png'],'Resolution',300)

end
