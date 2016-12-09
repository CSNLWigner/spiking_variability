function saveFigure(figureHandle,filename,removeLeadingZeros,panelLettering,groupings)    
    if nargin < 5
        groupings = {};
        if nargin < 4
            panelLettering = 'temporal';
            if nargin < 3
                removeLeadingZeros = true;
            end
        end
    end
    
    
    if ~strcmp(panelLettering,'none')
        if strcmp(panelLettering,'coords')
            letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';            
            xoffset = 0.12;
            xpanelwidth = 0.21;
            allfig = axes('position',[0 0 1 1],'Visible','off');
            firstrowY = 0.95;
            secondrowY = 0.48;
            ycoords = [firstrowY firstrowY firstrowY firstrowY firstrowY secondrowY secondrowY secondrowY secondrowY];
                        
            xcoords = [xoffset xoffset+xpanelwidth xoffset+2*xpanelwidth xoffset+3*xpanelwidth-0.04 0.84 xoffset xoffset+xpanelwidth xoffset+2*xpanelwidth xoffset+3*xpanelwidth];
            for l = 1:length(xcoords)
                text(xcoords(l), ycoords(l), letters(l), 'FontSize', 26, 'FontWeight', 'bold');
            end            
        else
            letterByPosition = strcmp(panelLettering,'position');
            letterSubplots(figureHandle,letterByPosition,groupings);
        end
    end
    if removeLeadingZeros
        removeLeadingZeroTicks(figureHandle);
    end
    
%     %allfig = axes('position',[0 0 1 1],'Visible','off');
%     text(xoffset+xpanelwidth-0.107, 0.108, 'H', 'FontSize', 30, 'FontWeight', 'bold', 'color','w');
%     text(xoffset+2*xpanelwidth-0.111, 0.108, 'H', 'FontSize', 30, 'FontWeight', 'bold', 'color','w');
    
    %path_name = ['figures/' filename];
    path_name = [filename];
    savefig(figureHandle,[path_name '.fig']);
    set(figureHandle,'PaperPositionMode','auto')
    print(figureHandle,[path_name '.png'],'-dpng','-r0');
    set(figureHandle,'Units','Inches');
    pos = get(figureHandle,'Position');
    set(figureHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(figureHandle,[path_name '.pdf'],'-dpdf','-r0');
end

function removeLeadingZeroTicks(figureHandle)
    figure(figureHandle);    
    childrenHandles = get(gcf,'children');
    for i = 1:length(childrenHandles)
        ah = childrenHandles(i);
        if ~strcmp(ah.Type,'axes')
            continue
        end
        nxtl = {};
        for j = 1:length(ah.XTickLabel)
           at =  ah.XTickLabel{j};
           if length(at) > 1
               if strcmp(at(1),'0')
                   at(1) = '';
               elseif strcmp(at(1),'-') && strcmp(at(2),'0')
                   at(2) = '';
               end
           end
           nxtl{end+1} = at;           
        end
        set(ah,'XTickLabel',nxtl)
        
        nytl = {};
        for j = 1:length(ah.YTickLabel)
           at =  ah.YTickLabel{j};
           if length(at) > 1
               if strcmp(at(1),'0')
                   at(1) = '';
               elseif strcmp(at(1),'-') && strcmp(at(2),'0')
                   at(2) = '';
               end
           end
           nytl{end+1} = at;           
        end
        set(ah,'YTickLabel',nytl)
    end
end
