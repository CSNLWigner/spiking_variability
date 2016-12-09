function equalizePlots(plotRows,plotColumns,minRow,maxRow,minColumn,maxColumn,equalizeX,equalizeY,relocateStars)
    if nargin < 9
        relocateStars = false;
    end
        
    minY = Inf;
    maxY = -Inf;
    minX = Inf;
    maxX = -Inf;
    for i = minRow:maxRow
        for j = minColumn:maxColumn
            subplot(plotRows,plotColumns,(i-1)*plotColumns+j);
            yl = ylim();
            xl = xlim();
            if yl(1) < minY
                minY = yl(1);
            end
            if yl(2) > maxY
                maxY = yl(2);
            end
            if xl(1) < minX
                minX = xl(1);
            end
            if xl(2) > maxX
                maxX = xl(2);
            end
        end
    end
    
    for i = minRow:maxRow
        for j = minColumn:maxColumn
            subplot(plotRows,plotColumns,(i-1)*plotColumns+j);
            if equalizeX
                xlim([minX maxX]);
            end
            if equalizeY
                ylim([minY maxY]);
                if relocateStars
                    h = findobj(gca,'Type','text');
                    for hi = 1:length(h)
                        h(hi).Position(2) = maxY - 0.07 * (maxY-minY);                    
                    end
                end
            end                        
        end
    end
end