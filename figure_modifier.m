clear all
close all
clc

%Figure settings
mkr_sz = 8;
line_wdth = 1.2;
fig_sz = [100 100 650 400];
manual_resize = 1;
tick_sz = [0.03 0.012];
font_sz = 20;
lgd = 1;
xlbl = 1;
ylbl = 1;
titl = 0;
xscl = 1;
yscl = 1;

sub_fig = 0;
n = 2;
m = 1;
naming_scheme = ["a)" "b)" "c)" "d)" "e)" "f)" "g)"];

make_gif = 1;

load("pumpFit_freq2velo.mat");

topfolder = '\figures\matlab_figs\';
if sub_fig == 1
    % [processing_files, processing_file_dir] = uigetfile(topfolder,'*.fig',sprintf('Select %d files'),n*m,'MultiSelect','on');
    subfig_fig = figure;
    subfig_tiled = tiledlayout(subfig_fig,n,m,'Padding','compact','TileSpacing','compact');
    num_figs = n*m;

else
    [temp, processing_file_dir] = uigetfile(topfolder,'*.fig','MultiSelect','on');
    % num_figs = length(temp);
    if class(temp) == 'char'
        processing_files = cellstr(temp);
        num_figs = 1;
    else
        processing_files = temp;
        num_figs = numel(processing_files);
    end
    
end

%%

for ii = 1:num_figs
    % close all
    if sub_fig == 1
        % ntile = nexttile;
        [processing_files, processing_file_dir] = uigetfile(topfolder,'*.fig',sprintf('Select fig #%d',ii),'MultiSelect','off');
        current_fig = openfig(fullfile(processing_file_dir,processing_files));
        figure(subfig_fig)
        subfig_subplot(ii) = subplot(n,m,ii);
        copyobj(allchild(get(current_fig,'CurrentAxes')),subfig_subplot(ii));
        % for ai = 1:numel(current_fig)
        %     newax = copyobj(current_fig(ai), subfig_tiled);
        %     % set(newax, 'Position', get(gca, 'Position')); % Align axes
        %     % delete(gca); % Remove placeholder axes
        % end
        lgd = input("Legend (0 no/1 yes)=");
        xlbl = input("x label (0 no/1 yes)=");
        ylbl = input("y label (0 no/1 yes)=");
        titl = input("title (0 no/1 yes)=");
        xscl = input("x scale (0 no/1 yes)=");
        yscl = input("y scale (0 no/1 yes)=");
    else
        current_fig = openfig(fullfile(processing_file_dir,processing_files{ii}));
    end

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
    set(gca,'OuterPosition',[0 0 0.9 1]);
    if ii==1 && manual_resize==1
        
        disp('Press key after figure has desired size')
        pause
        fig_sz = current_fig.Position;
        axes_sz = tightPosition(ax,Units="normalized");
    end
    
    current_fig.Position = fig_sz;
    ax.Position = axes_sz;

    if make_gif == 1
        h = flip(h);
        set(h,'Visible','off')
        for k = 1:numel(h)
            set(h(k),'Visible','on')
            drawnow
            exportgraphics(current_fig,['figures\slidefigs\' processing_files{ii}(1:end-4) '_' num2str(k) '.png'],'Resolution',400)
            
        end
    end

    exportgraphics(current_fig,['figures\' processing_files{ii}(1:end-4) '.pdf'],'Resolution',400)
    exportgraphics(current_fig,['figures\' processing_files{ii}(1:end-4) '.png'],'Resolution',400)

end
