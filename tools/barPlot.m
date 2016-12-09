function [ts,ps] = barPlot(data,labels,pvals,titleString,ylab,additional,errorquant,plotMedian,colors)    
    ylim([-0.04 0.07])
    xlim([0 4])
    plot([0 2.25],[0 0],'k')
    hold on
    plot([2.65 3.6],[0 0],'k')
    %plot([2.3 2.5],[-0.005 0.005],'k')
    %plot([2.5 2.7],[-0.005 0.005],'k')
    %allfig = axes('position',[0 0 1 1],'Visible','off');
    r = rectangle('Position',[2.25 -0.045 0.4 0.01],'FaceColor','w','EdgeColor','w');
    ylim([-0.04 0.07])
    set(r,'clipping','off')
    %text(-1,0,'H','FontSize',40)

    blue = [0 0.4470 0.7410];
    red = [1 0 0];
    grey = [0.7 0.7 0.7];

    means = arrayMean(data);
    medians = arrayMean(data,1,@median);
    sems = sem(data);
    nVars = length(data);

    if nargin < 9        
        colors = repmat(grey,nVars,1);        
        if nargin < 8
            plotMedian = false;
            if nargin < 7
                errorquant = 'SEM';
            end
        end
    end
    
    if isempty(colors)
        colors = repmat(red,nVars,1);        
    end
        
    % calculate p-values and assemble strings with them    
    pValStrings = '';
    nPs = size(pvals,1);
    
    CIs = zeros(nVars,2);
    for v1 = 1:nPs
        [~,act_p] = ttest2(data{pvals(v1,1)},data{pvals(v1,2)});
        labelParts1 = strsplit(labels{pvals(v1,1)});
        labelParts2 = strsplit(labels{pvals(v1,2)});
        pValStrings = [pValStrings sprintf(' p_{t-test}(%s,%s) = %.3f',labelParts1{1},labelParts2{1},act_p)];        
    end
    
    hold on
    trendLines = [];
    for v = 1:nVars
        if isempty(data{v})
            continue
        end
        if length(data{v}) > 1
            CIs(v,:) = dataConfinf(data{v},95);
        end

        if means(v)>medians(v)            
            %bar(v, means(v), 'facecolor', colors(v,:));
            if plotMedian
                bar(v, medians(v), 'facecolor', 1.2*colors(v,:));
            end
        else
            if plotMedian
                bar(v, medians(v), 'facecolor', 1.2*colors(v,:));
            end
            %bar(v, means(v), 'facecolor', colors(v,:));
        end        
        if length(data{v}) > 1
            plotConfinf(data{v},95,v,colors(v,:),false,false,'none',[means(v),sems(v)]);
        end
        if strcmp(additional,'scatter')
            scatter(repmat(v,length(data{v}),1) + randn(length(data{v}),1)/20,data{v},'filled','MarkerEdgeColor','k')                        
        elseif strcmp(additional,'trendLines')
            trendLines = [trendLines; data{v}'];        
        end

    end
    
    if strcmp(additional,'lineHist')
        ylims = ylim();
        for v = 1:nVars
            [actHist,binCenters] = hist(data{v},linspace(min(ylims(1),min(data{v})),max(ylims(2),max(data{v})),40));
            plot(smooth(repmat(v,length(actHist),1)+(actHist./(2*max(actHist)))',5),binCenters,'LineWidth',2)        
        end
        ylims(2) = max(ylims(2),(max(means)+max(sems))*1.2);
        ylim(ylims);
    elseif strcmp(additional,'trendLines')
        plot(trendLines);
    end
    
    if strcmp(errorquant,'SEM')
        error_to_plot = sems;
    elseif strcmp(errorquant,'CI')
        error_to_plot = CIs(:,1);
    end
    %errorbar(means,error_to_plot,'.','LineWidth',2,'Color',[0.1 0.1 0.1]);

    ylab
    ts = [];
    ps = [];
    if strcmp(additional,'zeroDiff')
        for v = 1:nVars
            if strcmp(errorquant,'CI')
                pval = min(sum(data{v} > 0),sum(data{v} < 0)) * (1 / length(data{v}));
                fprintf('%s p value: %.4f\n',labels{v},pval)
            else
                [~,pval,~,stats] = ttest(data{v});
                fprintf('%s p value: %.4f t stat: %.4f df: %d\n',labels{v},pval,stats.tstat,stats.df);
                ts = [ts;stats.tstat];
                ps = [ps;pval];
            end
            if pval < 0.05 && length(data{v}) > 1
                starstring = '*';
%                 if pval < 0.01
%                     starstring = '**';
%                 end
                yl = ylim();
                text(v, yl(2) - 0.07 * (yl(2)-yl(1)), starstring, 'FontSize', 30, 'HorizontalAlignment', 'center');                
            end
        end
    end
    
    xlabel(pValStrings);    
    if ~isempty(ylab)
        % ylabel(['$' ylab '$'],'Interpreter','LaTex');
        ylabel(ylab);
    end
    if ~isempty(titleString)
        title(titleString);
    end
    set(gca,'XTick',1:nVars,'XTickLabel',labels,'FontSize',16);
end

function m = arrayMean(data,dim,calcFunc)
    %params = parseVarargin(varargin,{'calcFunc'},{'mean'});
    if nargin < 3
        calcFunc = @mean;
    end
    
    if nargin < 2
        dim = 1;
    end
    if iscell(data)
        m = zeros(size(data));
        for i = 1:size(data,1)
            for j = 1:size(data,2)
                %m(i,j) = mean(data{i,j},dim);   
                if isempty(data{i,j})
                    m(i,j) = nan;
                else
                    m(i,j) = calcFunc(data{i,j},dim);
                end
            end
        end
    else
        %m = mean(data,dim);
        m = calcFunc(data,dim);
    end
end

function s = sem(data,dim)
    if nargin < 2
        dim = 1;
    end
    if iscell(data)
        s = zeros(size(data));
        for i = 1:size(data,1)
            for j = 1:size(data,2)
                if isempty(data{i,j})
                    s(i,j) = nan;
                else
                    s(i,j) = sem(data{i,j},dim);
                end
            end
        end
    else 
        s = std(data,0,dim) ./ sqrt(size(data,dim));
    end
end