function textHandles = letterSubplots(figureHandle,byPosition,groupings)
    if nargin < 3
        groupings = {};
        if nargin < 2
            byPosition = false;
        end
    end

    figure(figureHandle);    
    childrenHandles = get(gcf,'children');
    
    if byPosition
        orderedHandles = childrenHandles;
        for i = 1:length(orderedHandles)
            for j = 1:length(orderedHandles)-1
                if comparePositions(orderedHandles(j+1).Position, orderedHandles(j).Position)
                    temp = orderedHandles(j);
                    orderedHandles(j) = orderedHandles(j+1);
                    orderedHandles(j+1) = temp;
                end
            end
        end
    end
    
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    textHandles = [];
    for i = 1:length(childrenHandles)        
        if ~byPosition
            % this assumes subplots were added in order. otherwise we need to
            % order them based on positions
            ah = childrenHandles(end-i+1);
        else
            ah = orderedHandles(i);
        end
        if ~strcmp(ah.Type,'axes')
            continue
        end
        if ~strcmp(ah.Units,'normalized')
            error('Not implemented for position coordinates in units other than normalized');
        end
        axes(ah);
        % filter out insets
        if ah.Position(3) < 0.1 && ah.Position(4) < 0.15
            included = false;
%             for hi = 1:length(includeHandle)
%                 if ah == includeHandle{hi}
%                     included = true;
%                 end
%             end
            if ~included
                continue
            end
        end        
        withinGroup = false;
        firstInGroup = false;
        if ~isempty(groupings)
           for g = 1:length(groupings) 
               for e = 1:length(groupings{g})
                   if groupings{g}(e) == i
                       if e == 1
                           firstInGroup = true;
                       else
                           withinGroup = true;
                       end
                   end
               end
           end
        end
        
        xl = xlim();
        yl = ylim();        
        
        if withinGroup
            continue
        elseif firstInGroup
            % this is highly specific to one figure
            th = text(xl(1) - 0.8 * (xl(2)-xl(1)), yl(2) - 0.14 * (yl(2)-yl(1)), letters(1), 'FontSize', 26, 'FontWeight', 'bold');
        else
            th = text(xl(1) + 0.05 * (xl(2)-xl(1)), yl(2) + 0.07 * (yl(2)-yl(1)), letters(1), 'FontSize', 26, 'FontWeight', 'bold');
        end
        textHandles = [textHandles; th];
        letters(1) = '';
    end
end

function pos1first = comparePositions(pos1,pos2)
    if pos1(2) > pos2(2)
        pos1first = true; % it is more upwards
    elseif pos1(2) == pos2(2) && pos1(1) < pos2(1)
        pos1first = true; % it is more to the left
    else
        pos1first = false;
    end
end