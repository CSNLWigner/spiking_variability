function topDataAxis = plotRunHists(delta_std_ctr_P,delta_std_ctr_RG,varargin)
    parser = inputParser;
    addParameter(parser,'dataMinTrial',0,@isnumeric);
    addParameter(parser,'colors',[0.1 0.7]);
    addParameter(parser,'scattercol',[1 0 0]);
    addParameter(parser,'legend',false,@islogical);
    addParameter(parser,'xlabel',false,@islogical);
    addParameter(parser,'ylabel','',@ischar);
    parse(parser,varargin{:});
    params = parser.Results;
  
    if params.dataMinTrial > 0
        topDataAxis = plotDataDelta(params.dataMinTrial,4,4,2,1,params.scattercol)
        %plotDataDelta(params.dataMinTrial,1,1,1,1)
        %subplot('Position',[0.1 0.1 0.55 0.8])
        subplot(2,7,6)
        %subplot(1,2,1)
    end
    
    dsp_col = params.colors(1) * ones(1,3);
    rg_col = params.colors(2) * ones(1,3);
    data_limits = [-0.02 0.12];
    minsim = min([delta_std_ctr_P(:); delta_std_ctr_RG(:)]);
    maxsim = max([delta_std_ctr_P(:); delta_std_ctr_RG(:)]);
%     if params.dataMinTrial > 0
%         limits = [min([minsim data_limits(1) 0]) - 0.01 max([maxsim data_limits(2) 0]) + 0.01];
%     else
%         limits = [min([minsim 0]) - 0.01 max([maxsim 0]) + 0.01];
%     end
    
    limits = [min([minsim 0]) - 0.01 max([maxsim 0]) + 0.01];
    binedges = linspace(limits(1),limits(2),30);
    histogram(delta_std_ctr_P,binedges,'normalization','probability','FaceColor',dsp_col);
    hold on
    histogram(delta_std_ctr_RG,binedges,'normalization','probability','FaceColor',rg_col);
    yl_orig = ylim();
    plot([0 0],ylim() * 3,'k:','LineWidth',2)   
    ylim(yl_orig * 1.5);
    lgndStrs = {'DSP','RG'}; 
    box off
    
    %text(0,yl_orig(2) * 1.4,[sprintf('%d',sum(delta_std_ctr_P<0)) '\leftarrow'],'FontSize',16,'Color',dsp_col,'horizontalAlignment','right')
    %text(0,yl_orig(2) * 1.4,['\rightarrow' sprintf('%d',sum(delta_std_ctr_P>0))],'FontSize',16,'Color',dsp_col,'horizontalAlignment','left')
    %text(0,yl_orig(2) * 1.2,[sprintf('%d',sum(delta_std_ctr_RG<0)) '\leftarrow'],'FontSize',16,'Color',rg_col,'horizontalAlignment','right')
    %text(0,yl_orig(2) * 1.2,['\rightarrow' sprintf('%d',sum(delta_std_ctr_RG>0))],'FontSize',16,'Color',rg_col,'horizontalAlignment','left')
    
    plotConfinf(delta_std_ctr_P,0.95,yl_orig(2) * 1.1,dsp_col,true,true);
    plotConfinf(delta_std_ctr_RG,0.95,yl_orig(2) * 1.0,rg_col,true,true);
    
    if params.legend
        if abs(limits(1)) > abs(limits(2))
            lh = legend(lgndStrs,'Location','Northwest');
        else
            lh = legend(lgndStrs,'Location','Northeast');
        end
        legend boxoff
    end    
    xlim(limits)
    set(gca,'FontSize',16)    
    if params.xlabel
        xlabel('\Delta SD_{ctr}')
    end

    actPos = get(gca,'Position');
    set(gca,'Position',[0.735 actPos(2:4)]);
end