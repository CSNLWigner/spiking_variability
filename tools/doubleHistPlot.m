function doubleHistPlot(data1,data2,limits,resolution,legendStrings,figFilename,xlabels,colors)
    
    setColors = true;
    if nargin < 8
        setColors = false;
        colors = {'r','b'};
    end
    
    if nargin < 7
        xlabels = true;
    end

    if ~isempty(figFilename)
        figHist = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
    end        
    
    if ~isempty(limits) 
        minVal = limits(1);
        maxVal = limits(2);
    end

    if isempty(data2)
        if isempty(limits) 
            minVal = min(data1);
            maxVal = max(data1);
        end
        hist(data1,linspace(minVal,maxVal,resolution));
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor','k','EdgeColor','w','facealpha',0.35)
    else
        if isempty(limits) 
            minVal = min(min(data1),min(data2));
            maxVal = max(max(data1),max(data2));
        end
        hi = histogram(data1,linspace(minVal,maxVal,resolution));        
        hi.Normalization = 'probability';        
        hold on
        hi2 = histogram(data2,linspace(minVal,maxVal,resolution));
        hi2.Normalization = 'probability';
        if setColors
            hi.FaceColor = colors{1};
            hi2.FaceColor = colors{2};
        end
%         hi.FaceAlpha = 0.3;
%         hi2.FaceAlpha = 0.3;
       
%         %hist(data1,linspace(minVal,maxVal,resolution));
%         %h = findobj(gca,'Type','patch');
%         %set(h,'FaceColor',colors{1},'EdgeColor','w','facealpha',0.75)
%         %set(h,'facealpha',0.75)
%         hold on
%         hist(data2,linspace(minVal,maxVal,resolution));
%         hold off
%         h = findobj(gca,'Type','patch');
%         %set(h,'facealpha',0.75);
%         %set(h,'FaceColor',colors{2},'EdgeColor','w','facealpha',0.75)
%         set(h,'EdgeColor','w','facealpha',0.75)
        if ~isempty(legendStrings)
            legend(legendStrings);
        end
    end
        
    xlim([minVal maxVal])
    ylim([0 0.2])
    set(gca,'Ytick',[],'FontSize',16, 'box','off');
    
    [histVals,binCenters] = hist(data1,linspace(minVal,maxVal,resolution));
    [~,maxIdx] = max(histVals);
    mode1 = binCenters(maxIdx);
    
    [histVals,binCenters] = hist(data2,linspace(minVal,maxVal,resolution));
    [~,maxIdx] = max(histVals);
    mode2 = binCenters(maxIdx);
    
    if xlabels
        xlabel({['\color{red}'  sprintf('mean %.3f median %.3f mode %.3f std %.3f krt %.3f',mean(data1),median(data1),mode1,std(data1),kurtosis(data1))], ...
                ['\color{blue}' sprintf('mean %.3f median %.3f mode %.3f std %.3f krt %.3f',mean(data2),median(data2),mode2,std(data2),kurtosis(data2))]});
                %['\color{black}' xlabelString]});
    end

    if ~isempty(figFilename)
        set(gcf,'PaperPositionMode','auto')    
        print(figHist,figFilename,'-dpng','-r0');
        close(figHist);
    end
end