function gen_fig_9(recalc)
    nPriorSample = 9;
    nUnit = 50;
    corrbeta_P = 1.5;
    corrbeta_RG = 6;
    nBin = 1;
    nRerun = 100;
    nTrial = 1000;    
    %contrast_mean_factor = 2.1;
    contrast_var_factor = 0.7;
    
    
    nHistPlots = 6;
    histPlotIndices = floor(linspace(1,nPriorSample,nHistPlots));
    grey_p = 0.1;
    grey_r = 0.7;
    
%     meta_tc_P = [50 30 30 70 70];
%     meta_tc_RG = [80 60 60 100 100];
%     meta_var_P = [4 2 6 2 6];
%     meta_var_RG = [4 2 6 2 6];
    
    meta_mean_ctr_dep = [1.5 2 2.5];
%     meta_tc_P = 50 * ones(1,length(meta_mean_ctr_dep));
%     meta_tc_RG = 80 * ones(1,length(meta_mean_ctr_dep));
    meta_tc_P = 60 * ones(1,length(meta_mean_ctr_dep));
    meta_tc_RG = 130 * ones(1,length(meta_mean_ctr_dep));
    meta_var_P = 4 * ones(1,length(meta_mean_ctr_dep));
    meta_var_RG = 4 * ones(1,length(meta_mean_ctr_dep));
    
    descr = {'','Lower mean, lower var in MP','Lower mean, higher var in MP','Highrt mean, lower var in MP','Higher mean, higher var in MP'};
    
    nMeta = length(meta_tc_P);
    nCol = 1;
    plotshift = 0.01;
    close all
    
%     [all_spikes, conditions] = loadEckerData(38);
%     [ecker_sc_lc, ecker_sc_hc] = separateContrast(all_spikes, conditions);    
%     [cEsc_lc,~,~] = getStats(ecker_sc_lc);
%     [cEsc_hc,~,~] = getStats(ecker_sc_hc);          
%     delta_std_ctr_E = zeros(nRerun,1);
%     cEsc_hc_sh = cEsc_hc(randperm(length(cEsc_hc)));
%     cEsc_lc_sh = cEsc_lc(randperm(length(cEsc_lc)));
%     bootsize_hc = floor(length(cEsc_hc)/nRerun);
%     bootsize_lc = floor(length(cEsc_lc)/nRerun);
%      for b = 1:nRerun            
%         endindex_hc = b*bootsize_hc;
%         endindex_lc = b*bootsize_lc;
%         if b==nRerun
%             endindex_hc = length(cEsc_hc);
%             endindex_lc = length(cEsc_lc);
%         end
%         cEsc_hc_act = cEsc_hc_sh((b-1)*bootsize_hc + 1 : endindex_hc);
%         cEsc_lc_act = cEsc_lc_sh((b-1)*bootsize_lc + 1 : endindex_lc);
%         delta_std_ctr_E(b) = std(cEsc_hc_act)-std(cEsc_lc_act);
%      end
    
    if recalc
        allStats = cell(nMeta,4);
        for m = 1:nMeta
            fprintf('Metaparams %d/%d\n ',m,nMeta);

        %     allPvals = zeros(nPriorSample,2);
        %     allDeltaMeans = zeros(nPriorSample,2);    
        %     dataDeltas = zeros(nPriorSample,2);
        %     dataPvals = zeros(nPriorSample,2);    
            allCorrDeltas = zeros(nPriorSample,nRerun,2);
            allRates = zeros(nPriorSample,nRerun,2,2);
            allFFs = zeros(nPriorSample,nRerun,2,2);

            %parameterSamples = rand(nPriorSample,1) + 1; 
            parameterSamples = linspace(0.8,1.2,nPriorSample); 

            for s = 1:nPriorSample
                if nPriorSample >= nRerun
                    printCounter(s, 'maxVal', nPriorSample, 'stringVal', 'Prior sample');
                else
                    fprintf('Prior sample %d/%d ',nPriorSample,s);
                end
                delta_std_ctr_P = zeros(nRerun,1);
                delta_std_ctr_RG = zeros(nRerun,1);

                lc_mean_rate_P = zeros(nRerun,1);
                lc_mean_rate_RG = zeros(nRerun,1);
                lc_mean_ff_P = zeros(nRerun,1);
                lc_mean_ff_RG = zeros(nRerun,1);

                hc_mean_rate_P = zeros(nRerun,1);
                hc_mean_rate_RG = zeros(nRerun,1);
                hc_mean_ff_P = zeros(nRerun,1);
                hc_mean_ff_RG = zeros(nRerun,1);

                for i = 1:nRerun
                    if nPriorSample < nRerun
                        printCounter(i, 'maxVal', nRerun, 'stringVal', 'Rerun');
                    end
                    if s == 1 && i == 1
                        corrparam_P = corrbeta_P;
                        corrparam_RG = corrbeta_RG;
                    else
        %                 corrparam_P = corrmat_P;
        %                 corrparam_RG = corrmat_RG;
                    end
                    %[cPsc_lc_act,cPsc_hc_act,cRGsc_lc_act,cRGsc_hc_act,~,~,~,~,corrmat_P,corrmat_RG,~] = simulateCorrelations(nUnit, corrparam_P, corrparam_RG, nBin,'rate_exponent',rateExponents(s),'n_trial',nTrial);
                    %[Pstat,RGstat] = simulateCorrelations(nUnit, corrparam_P, corrparam_RG, nBin,'rate_exponent',rateExponents(s),'n_trial',nTrial);
                    contrast_corr_factor = parameterSamples(s);
                    [Pstat,RGstat] = simulateCorrelations(nUnit, corrparam_P, corrparam_RG, nBin, ...
                        'rate_exponent',1.4,'n_trial',nTrial,'contrast_corr_factor',contrast_corr_factor,'contrast_mean_factor',meta_mean_ctr_dep(m),'contrast_var_factor',contrast_var_factor, ...
                        'mp_var_b_P',meta_var_P(m),'tc_coeff_P',meta_tc_P(m),'mp_var_b_RG',meta_var_RG(m),'tc_coeff_RG',meta_tc_RG(m));
        %             corrmat_P = Pstat.HC.MP.corrmat;
        %             corrmat_RG = RGstat.HC.MP.corrmat;
                    % 'contrast_corr_factor_P',1,'contrast_corr_factor_RG',1,'mp_exponent',1,            
                    delta_std_ctr_P(i) = std(Pstat.HC.SC.corrvec)-std(Pstat.LC.SC.corrvec);
                    delta_std_ctr_RG(i) = std(RGstat.HC.SC.corrvec)-std(RGstat.LC.SC.corrvec);

                    lc_mean_rate_P(i) = mean(Pstat.LC.SC.meanvec);
                    lc_mean_rate_RG(i) = mean(RGstat.LC.SC.meanvec);
                    lc_mean_ff_P(i) = mean(Pstat.LC.SC.ffvec);
                    lc_mean_ff_RG(i) = mean(RGstat.LC.SC.ffvec);

                    hc_mean_rate_P(i) = mean(Pstat.HC.SC.meanvec);
                    hc_mean_rate_RG(i) = mean(RGstat.HC.SC.meanvec);
                    hc_mean_ff_P(i) = mean(Pstat.HC.SC.ffvec);
                    hc_mean_ff_RG(i) = mean(RGstat.HC.SC.ffvec);

        %             delta_std_ctr_P(i) = std(cPsc_hc_act)-std(cPsc_lc_act);
        %             delta_std_ctr_RG(i) = std(cRGsc_hc_act)-std(cRGsc_lc_act);
                end        

        %         allDeltaMeans(s,1) = mean(delta_std_ctr_P);
        %         allDeltaMeans(s,2) = mean(delta_std_ctr_RG);

                allCorrDeltas(s,:,1) = delta_std_ctr_P;
                allCorrDeltas(s,:,2) = delta_std_ctr_RG;
                allRates(s,:,1,1) = lc_mean_rate_P;
                allRates(s,:,2,1) = lc_mean_rate_RG;
                allRates(s,:,1,2) = hc_mean_rate_P;
                allRates(s,:,2,2) = hc_mean_rate_RG;

                allFFs(s,:,1,1) = lc_mean_ff_P;
                allFFs(s,:,2,1) = lc_mean_ff_RG;
                allFFs(s,:,1,2) = hc_mean_ff_P;
                allFFs(s,:,2,2) = hc_mean_ff_RG;

        %         dataDeltas(s,1) = mean(delta_std_ctr_P) - mean(delta_std_ctr_E);
        %         dataDeltas(s,2) = mean(delta_std_ctr_RG) - mean(delta_std_ctr_E);

        %         [~,allPvals(s,1),~,~] = ttest(delta_std_ctr_P);
        %         [~,allPvals(s,2),~,~] = ttest(delta_std_ctr_RG);
        %         
        %         [~,dataPvals(s,1),~,~] = ttest2(delta_std_ctr_P,delta_std_ctr_E);
        %         [~,dataPvals(s,2),~,~] = ttest2(delta_std_ctr_RG,delta_std_ctr_E);



                %allPvals(s,:) = poisson_integ_population(100,1.5,8,'makePlots',false,'bulkSim',true,'simVariedMPCorr',false,'rate_exponent',rand()+1,'contrast_corr_factor',1);
            end
            
            allStats{m,1} = allCorrDeltas;
            allStats{m,2} = allRates;
            allStats{m,3} = allFFs;
            allStats{m,4} = parameterSamples;
        end
        save('saveStats.mat','allStats')
    else
        load('saveStats.mat')
    end
    
    figCorrFact = figure('units','normalized','outerposition',[0 0 1 1]);
    for m = 1:nMeta
        
        figure(figCorrFact);

        %labelstring = 'Rate exponent';
        labelstring = 'Contrast-dependence of corr.';
    %     
    %     subplot(nRow,2,1)
    %     histogram(allDeltaMeans(:,1),linspace(-0.02,0.02,50),'normalization','probability')
    %     hold on
    %     histogram(allDeltaMeans(:,2),linspace(-0.02,0.02,50),'normalization','probability')
    %     hold off
    %     xlabel(['Mean of \Delta SD over ' sprintf('%d reruns',nRerun)]);
    %     
    %     subplot(nRow,2,2)
    %     histogram(allPvals(:,1),logspace(-5,0,50),'normalization','probability')
    %     hold on
    %     histogram(allPvals(:,2),logspace(-5,0,50),'normalization','probability')
    %     hold off
    %     xlabel('Nonzero \Delta SD P-value')
    %     set(gca,'XScale','log');
    %     legend({'P','RG'})
        parameterSamples = allStats{m,4};
        %subplot(nMeta,nCol,(m-1)*nCol+1)
        subplot(nMeta,nCol,(nMeta -m)*nCol+1)
        scatter(repmat(parameterSamples,[1 nRerun]),reshape(allStats{m,1}(:,:,1),[nPriorSample*nRerun 1]),[],grey_p*ones(nPriorSample*nRerun,3))
        hold on        
        scatter(repmat(parameterSamples,[1 nRerun]) + plotshift,reshape(allStats{m,1}(:,:,2),[nPriorSample*nRerun 1]),[],grey_r*ones(nPriorSample*nRerun,3))
        yl = ylim();
        hold off
        plot([min(parameterSamples)-plotshift max(parameterSamples)+plotshift],[0 0],'k--')
        
        ylim(yl);
        for s = 1:nPriorSample
            P_deltas = allStats{m,1}(s,:,1);
            RG_deltas = allStats{m,1}(s,:,2);
            [~,h_p] = plotConfinf(P_deltas,0.95,parameterSamples(s)-plotshift/2,grey_p*ones(1,3),false,false,'pval',[],true);
            [~,h_rg] = plotConfinf(RG_deltas,0.95,parameterSamples(s)+plotshift/2,grey_r*ones(1,3),false,false,'pval',[],true);
        end
        
        xlim([min(parameterSamples)-plotshift max(parameterSamples)+plotshift])        
        hold off
        
        if m==1            
            %title(labelstring);
            %legend([h_p h_rg],{'DSP','RG'},'FontSize',16,'Location','Northwest')
            xlabel('correlation factor')
        elseif m == nMeta
            %xlabel('correlation factor')
            legend([h_p h_rg],{'DSP','RG'},'FontSize',16,'Location','Northwest')
        end
        set(gca,'FontSize',16)
        %ylabel(sprintf('mean factor = %.2f',meta_mean_ctr_dep(m)));        
        ylabel('\Delta SD_{ctr}')
        
%         subplot(nMeta,nCol,(m-1)*nCol+2)        
%         scatter(repmat(parameterSamples,[1 nRerun]),reshape(allStats{m,2}(:,:,1,1),[nPriorSample*nRerun 1]))
%         hold on
%         scatter(repmat(parameterSamples,[1 nRerun]) + plotshift,reshape(allStats{m,2}(:,:,1,2),[nPriorSample*nRerun 1]))
%         scatter(repmat(parameterSamples,[1 nRerun]) + 2*plotshift,reshape(allStats{m,2}(:,:,2,1),[nPriorSample*nRerun 1]))
%         scatter(repmat(parameterSamples,[1 nRerun]) + 3*plotshift,reshape(allStats{m,2}(:,:,2,2),[nPriorSample*nRerun 1]))
%         %plot([min(parameterSamples) max(parameterSamples)+plotshift],[0 0],'k--')
%         xlim([min(parameterSamples) max(parameterSamples)+3*plotshift])
%         hold off
%         if m==1
%             legend({'P LC','P HC','RG LC','RG HC'})            
%         end
%         xlabel(sprintf('TC P %d TC RG %d V P %d V RG %d',meta_tc_P(m),meta_tc_RG(m),meta_var_P(m),meta_var_RG(m)))
%         
%         subplot(nMeta,nCol,(m-1)*nCol+3)        
%         scatter(repmat(parameterSamples,[1 nRerun]),reshape(allStats{m,3}(:,:,1,1),[nPriorSample*nRerun 1]))
%         hold on
%         scatter(repmat(parameterSamples,[1 nRerun]) + plotshift,reshape(allStats{m,3}(:,:,1,2),[nPriorSample*nRerun 1]))
%         scatter(repmat(parameterSamples,[1 nRerun]) + 2*plotshift,reshape(allStats{m,3}(:,:,2,1),[nPriorSample*nRerun 1]))
%         scatter(repmat(parameterSamples,[1 nRerun]) + 3*plotshift,reshape(allStats{m,3}(:,:,2,2),[nPriorSample*nRerun 1]))
%         %plot([min(parameterSamples) max(parameterSamples)+plotshift],[0 0],'k--')
%         xlim([min(parameterSamples) max(parameterSamples)+3*plotshift])
%         hold off
%         if m==1
%             legend({'P LC','P HC','RG LC','RG HC'})
%         end
        

    %     subplot(nRow,2,2)
    %     scatter(parameterSamples,allPvals(:,1))
    %     hold on
    %     scatter(parameterSamples,allPvals(:,2))
    %     hold off
    %     xlabel(labelstring);
    %     ylabel('Nonzero \Delta SD P-value')
    %     set(gca,'YScale','log');

    %     subplot(nRow,2,5)
    %     scatter(parameterSamples,dataDeltas(:,1))
    %     hold on
    %     scatter(parameterSamples,dataDeltas(:,2))
    %     plot([1 2],[0 0],'k--')
    %     hold off
    %     xlabel(labelstring);
    %     ylabel(['Diff of \Delta SD from data ' sprintf('%d reruns',nRerun)]);
    %     
    %     subplot(nRow,2,6)
    %     scatter(parameterSamples,dataPvals(:,1))
    %     hold on
    %     scatter(parameterSamples,dataPvals(:,2))
    %     hold off
    %     xlabel(labelstring);
    %     ylabel('Different from data \Delta SD P-value')
    %     set(gca,'YScale','log');                
    
    end
    
    saveFigure(figCorrFact,'fig9_corr_factors',false,'none');                                                        
    
%     figure('units','normalized','outerposition',[0 0 1 0.5]);
%     for m = 1:nMeta
%         ylab = sprintf('mean factor %.2f',meta_mean_ctr_dep(m));
%         for h = 1:nHistPlots
%             if h > 1
%                 ylab = '';
%             end
%             subplot(nMeta,nHistPlots,(m-1)*nHistPlots+h)
%             plotRunHists(allStats{m,1}(histPlotIndices(h),:,1),allStats{m,1}(histPlotIndices(h),:,2), ...
%                 'colors',[grey_p grey_r],'legend',m==1&&h==1,'xlabel',m==nMeta&&h==ceil(nHistPlots/2),'ylabel',ylab)
%             if h==1 && m==1
%                 title(sprintf('corr factor %.2f',parameterSamples(histPlotIndices(h))));
%             else
%                 title(sprintf('%.2f',parameterSamples(histPlotIndices(h))));
%             end
%             xl = round(xlim*100)/100;
%             xlim(xl);
%             set(gca,'XTick',[xl(1) 0 xl(2)]);
% 
%             if h > 1
%                 set(gca,'YTick',[]);
%             end
%         end
%     end
%     equalizePlots(nMeta,nHistPlots,1,nMeta,1,nHistPlots,false,true,false);
end
