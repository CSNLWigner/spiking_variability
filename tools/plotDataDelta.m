function firstAxisHandle = plotDataDelta(minTrials,nHorPanels,horPanelNo,nVerPanels,verPanelNo,col)
    load('ecker_data_v1_binned_static.mat')
    [all_spikes, conditions, trialNum, contrastLevels, unitNums] = loadEckerData(minTrials);
    pairNums = unitNums .* (unitNums - 1) ./ 2;
    nSess = length(all_spikes);
    hc_sd = zeros(nSess,1);
    lc_sd = zeros(nSess,1);    
    sd_ctr_trialnum = zeros(nSess*2,4);
                
    for i = 1:nSess
        [ecker_sc_lc, ecker_sc_hc] = separateContrast(all_spikes(i), conditions); 
        [cEsc_lc,~,~] = getStats(ecker_sc_lc);
        [cEsc_hc,~,~] = getStats(ecker_sc_hc); 
        hc_sd(i) = std(cEsc_hc);
        lc_sd(i) = std(cEsc_lc);
        sd_ctr_trialnum(i,:) = [std(cEsc_hc) size(data{i}.spikes,4) data{i}.conditions(2).contrast mean(cEsc_hc)];
        sd_ctr_trialnum(nSess+i,:) = [std(cEsc_lc) size(data{i}.spikes,4) data{i}.conditions(1).contrast mean(cEsc_lc)];
    end
    %sd_ctr_trialnum
%     subplot(1,2,1)
%     scatter(sd_ctr_trialnum(:,3),sd_ctr_trialnum(:,1),sd_ctr_trialnum(:,2).^2)
%     subplot(1,2,2)
%     scatter(sd_ctr_trialnum(:,3),sd_ctr_trialnum(:,4),sd_ctr_trialnum(:,2).^2)
    
%     ctrdiffIntensity = sqrt(1-((contrastLevels(:,2)-contrastLevels(:,1)) / 100));
%     %colors = repmat(trialIntensity, [1 3]);
%     colors = [0.7 * ctrdiffIntensity 0.9 * ctrdiffIntensity 0.7 * ctrdiffIntensity];    
%     unitNumIntensity = 1 - pairNums ./ max(pairNums);
%     unitColors = [1-unitNumIntensity unitNumIntensity unitNumIntensity];
%     scatter(hc_sd - lc_sd, 0.05*ones(nSess,1), trialNum.^3 / 500, colors, 'filled','MarkerEdgeColor','k','LineWidth',2);
%     hold on
%     scatter(hc_sd - lc_sd, 0.05*ones(nSess,1), trialNum.^3 / 500, unitColors,'LineWidth',2);
    
    fs = 12;
    left = (0.75 / nHorPanels) + (horPanelNo-1)/nHorPanels;
    width = 0.18 / nHorPanels;
    verPanelNo = nVerPanels - verPanelNo + 1;
    height = 0.18 / nVerPanels;
    vspace = 0.15 / nVerPanels;
    vstart = (verPanelNo-1)/nVerPanels;
    heightCorrection = 0.8;
    
    %subplot('Position',[left vstart+vspace width height])
    firstAxisHandle = subplot(7,12,12);
    %subplot(3,3,3);
    scatter(contrastLevels(:,2)-contrastLevels(:,1),hc_sd - lc_sd,'MarkerFaceColor',col,'MarkerEdgeColor',col)
    xlabel('Ctr diff')
    set(gca,'FontSize',fs)
    actPos = get(gca,'Position');set(gca,'Position',[actPos(1:3) actPos(4)*heightCorrection]);
    
    %subplot('Position',[left vstart+2*vspace+height width height])
    subplot(7,12,24);
    %subplot(3,3,6);
    scatter(trialNum,hc_sd - lc_sd,'MarkerFaceColor',col,'MarkerEdgeColor',col)
    xlabel('No. of trials')
    yLab = ylabel('Experiment \Delta SD_{ctr}');
    yLab.Position(1) = yLab.Position(1) + 5;
    set(gca,'FontSize',fs)
    actPos = get(gca,'Position');set(gca,'Position',[actPos(1:3) actPos(4)*heightCorrection]);
    
    %subplot('Position',[left vstart+3*vspace+2*height width height])
    subplot(7,12,36);
    %subplot(3,3,9);
    scatter(pairNums,hc_sd - lc_sd,'MarkerFaceColor',col,'MarkerEdgeColor',col)
    xlabel('No. of pairs')
    set(gca,'FontSize',fs)
    actPos = get(gca,'Position');set(gca,'Position',[actPos(1:3) actPos(4)*heightCorrection]);
end
