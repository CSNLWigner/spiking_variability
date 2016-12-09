function [CI,h,pval] = plotConfinf(x,conflevel,other_coord,color,horizontal,keepLimits,sig,precomp_mu_CI,wideFigure)
    if nargin < 9
        wideFigure = false;
    end
    if nargin < 8 || isempty(precomp_mu_CI)
        [m,CI] = dataConfinf(x,conflevel);
    else
        m = precomp_mu_CI(1);
        CI = [m-precomp_mu_CI(2) m+precomp_mu_CI(2)];
    end
    
%     [m,~,CI,~] = normfit(x,conflevel);
    
    yl = ylim();
    xl = xlim();
    yrange = yl(2) - yl(1);
    xrange = xl(2) - xl(1);
    hold on
    pval = Inf;
    
    if horizontal
        h = plot(CI,other_coord*ones(2,1),'Color',color,'LineWidth',2);
        meantick = yrange*0.02;
        plot(m*ones(2,1),[other_coord - meantick ,other_coord + meantick],'Color',color,'LineWidth',2)
        endtick = meantick;
        plot(CI(1)*ones(2,1),[other_coord - endtick ,other_coord + endtick],'Color',color,'LineWidth',2)
        plot(CI(2)*ones(2,1),[other_coord - endtick ,other_coord + endtick],'Color',color,'LineWidth',2)
        brackettick = xrange * 0.01;
        plot([CI(1) CI(1)+brackettick],(other_coord + endtick) * ones(2,1),'Color',color,'LineWidth',2)
        plot([CI(1) CI(1)+brackettick],(other_coord - endtick) * ones(2,1),'Color',color,'LineWidth',2)
        plot([CI(2)-brackettick CI(2)],(other_coord + endtick) * ones(2,1),'Color',color,'LineWidth',2)
        plot([CI(2)-brackettick CI(2)],(other_coord - endtick) * ones(2,1),'Color',color,'LineWidth',2)
    else
        h = plot(other_coord*ones(2,1),CI,'Color',color,'LineWidth',2);   
        
        if wideFigure
            meantick = xrange*0.003; % this is for the screen-wide figures
        else
            meantick = xrange*0.03; % this is for the 3-bar figures
        end
        plot([other_coord - meantick ,other_coord + meantick],m*ones(2,1),'Color',color,'LineWidth',2)
        endtick = meantick;
        plot([other_coord - endtick ,other_coord + endtick],CI(1)*ones(2,1),'Color',color,'LineWidth',2)
        plot([other_coord - endtick ,other_coord + endtick],CI(2)*ones(2,1),'Color',color,'LineWidth',2)
        brackettick = yrange * 0.005;
        plot((other_coord + endtick) * ones(2,1),[CI(1) CI(1)+brackettick],'Color',color,'LineWidth',2)
        plot((other_coord - endtick) * ones(2,1),[CI(1) CI(1)+brackettick],'Color',color,'LineWidth',2)
        plot((other_coord + endtick) * ones(2,1),[CI(2)-brackettick CI(2)],'Color',color,'LineWidth',2)
        plot((other_coord - endtick) * ones(2,1),[CI(2)-brackettick CI(2)],'Color',color,'LineWidth',2)
        
        %[~,pval] = ttest(x);
        textdist = 5*brackettick;
        if m < 0
            textpos = CI(1) - textdist;
            pval = sum(x > 0) * (1 / length(x));
            marker = 'd';
        else
            textpos = CI(2) + textdist;
            pval = sum(x < 0) * (1 / length(x));
            marker = 'o';
        end
        if nargin < 7 || strcmp(sig,'pval')
            text(other_coord,textpos,sprintf('%.2f',pval),'FontSize',12,'Color',color,'horizontalAlignment','center')
        elseif strcmp(sig,'star')
            % TODO do smthg if horizontal
            if pval < 0.05   
                text(other_coord, yl(2) - 0.07 * (yl(2)-yl(1)), '*', 'FontSize', 30, 'HorizontalAlignment', 'center');
            end
        end
        markersize = 100;
        if pval < 0.05   
            markercolor = color;
        else
            markercolor = 'w';
        end
        %scatter(other_coord,yl(1)+textdist,markersize,marker,'MarkerFaceColor',markercolor,'MarkerEdgeColor','k');        
    end
    if nargin < 6 || keepLimits
        ylim(yl)
        xlim(xl)
    end
end