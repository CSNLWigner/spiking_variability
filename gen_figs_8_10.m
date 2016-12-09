function pVals = gen_figs_8_10(n_unit, corrbeta_P, corrbeta_RG, varargin)
        
    parser = inputParser;
    addParameter(parser,'seed','shuffle');
    addParameter(parser,'makePlots',true,@islogical);
    addParameter(parser,'bulkSim',false,@islogical);
    addParameter(parser,'loadSavedSim',false,@islogical);
    addParameter(parser,'simVariedMPCorr',true,@islogical);
    
    addParameter(parser,'nRerun',100,@isnumeric);
    addParameter(parser,'nRerunVar',5,@isnumeric);
    addParameter(parser,'nBootstrap',5,@isnumeric);
    addParameter(parser,'nBinLong',20,@isnumeric);
    addParameter(parser,'nBinShort',1,@isnumeric);
    addParameter(parser,'minDataTrials',38,@isnumeric);
    addParameter(parser,'DSPBetaRange',[0.5 1 2 5 10]);
    addParameter(parser,'nMPVarSteps',3,@isnumeric);
    
    addParameter(parser,'contrast_corr_factor',1,@isnumeric);    
    addParameter(parser,'contrast_mean_factor',2,@isnumeric);    
    addParameter(parser,'contrast_var_factor',0.7,@isnumeric);    
    addParameter(parser,'mp_exponent',1,@isnumeric);
    addParameter(parser,'rate_exponent',1.4,@isnumeric);
    addParameter(parser,'tc_coeff_P',60,@isnumeric);    
    addParameter(parser,'tc_coeff_RG',130,@isnumeric);
    addParameter(parser,'mp_var_b_P',4,@isnumeric);
    addParameter(parser,'mp_var_b_RG',4,@isnumeric);    
    
    parse(parser,varargin{:});        
    params = parser.Results;

    setrandseed(params.seed);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [all_spikes, conditions, ~,~,~, dataBinWidthSec] = loadEckerData(params.minDataTrials);

    data_unit_num = 0;
    for i = 1:length(all_spikes)
        data_unit_num = data_unit_num + size(all_spikes{i},1);
    end
            
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    simCorr = @(nUnit,corrparam_P,corrparam_RG,nBin,baselineMean_P,baselineVariance_P) ...
        simulateCorrelations(n_unit, corrparam_P, corrparam_RG, nBin, 'tc_coeff_P',baselineMean_P, 'mp_var_b_P',baselineVariance_P, ...
        'contrast_corr_factor',params.contrast_corr_factor,'contrast_mean_factor',params.contrast_mean_factor,'contrast_var_factor',params.contrast_var_factor, ...
        'mp_exponent',params.mp_exponent,'rate_exponent',params.rate_exponent,'tc_coeff_RG',params.tc_coeff_RG);
    
    if ~params.loadSavedSim
        [Pstat,RGstat,simSmallBin] = simCorr(n_unit, corrbeta_P, corrbeta_RG, params.nBinLong, params.tc_coeff_P, params.mp_var_b_P);
        delta_std_ctr_P = zeros(params.nRerun,1);
        delta_std_ctr_RG = zeros(params.nRerun,1);
        delta_std_ori_P = zeros(params.nRerun,1);
        delta_std_ori_RG = zeros(params.nRerun,1);
        delta_mean_ctr_P = zeros(params.nRerun,1);
        delta_mean_ctr_RG = zeros(params.nRerun,1);
        delta_mean_ori_P = zeros(params.nRerun,1);
        delta_mean_ori_RG = zeros(params.nRerun,1);
        for i = 1:params.nRerun
            printCounter(i, 'maxVal', params.nRerun, 'stringVal', 'Errorbar');            
            [Pstat_act,RGstat_act] = simCorr(data_unit_num, corrbeta_P, corrbeta_RG, params.nBinShort, params.tc_coeff_P, params.mp_var_b_P);
            delta_std_ctr_P(i) = std(Pstat_act.HC.SC.corrvec)-std(Pstat_act.LC.SC.corrvec);
            delta_std_ctr_RG(i) = std(RGstat_act.HC.SC.corrvec)-std(RGstat_act.LC.SC.corrvec);
            delta_std_ori_P(i) = std(Pstat_act.HC.PO.SC.corrvec)-std(Pstat_act.HC.NO.SC.corrvec);
            delta_std_ori_RG(i) = std(RGstat_act.HC.PO.SC.corrvec)-std(RGstat_act.HC.NO.SC.corrvec);

            delta_mean_ctr_P(i) = mean(Pstat_act.HC.SC.corrvec)-mean(Pstat_act.LC.SC.corrvec);
            delta_mean_ctr_RG(i) = mean(RGstat_act.HC.SC.corrvec)-mean(RGstat_act.LC.SC.corrvec);
            delta_mean_ori_P(i) = mean(Pstat_act.HC.PO.SC.corrvec)-mean(Pstat_act.HC.NO.SC.corrvec);
            delta_mean_ori_RG(i) = mean(RGstat_act.HC.PO.SC.corrvec)-mean(RGstat_act.HC.NO.SC.corrvec);
        end
        
        cbp = params.DSPBetaRange;
        mp_std = zeros(length(cbp), params.nRerunVar);
        ori_dstd = zeros(length(cbp), params.nRerunVar);
        ctr_dstd = zeros(length(cbp), params.nRerunVar);
        
        varianceSteps = params.nMPVarSteps;
        ff_p = zeros(varianceSteps, params.nRerunVar);
        ori_dstd_varstep = zeros(varianceSteps, params.nRerunVar);
        ctr_dstd_varstep = zeros(varianceSteps, params.nRerunVar);
        
        if params.simVariedMPCorr
            for b = 1:length(cbp)        
                for r = 1:params.nRerunVar
                    printCounter(r, 'maxVal', params.nRerunVar, 'stringVal', sprintf('Beta %.2f',cbp(b)));
                    [Pstat_act,~] = simCorr(data_unit_num, cbp(b), 0, params.nBinShort, params.tc_coeff_P, params.mp_var_b_P);
                    mp_std(b,r) = std(Pstat_act.HC.MP.corrvec);
                    ori_dstd(b,r) = std(Pstat_act.HC.PO.SC.corrvec) - std(Pstat_act.HC.NO.SC.corrvec);
                    ctr_dstd(b,r) = std(Pstat_act.HC.SC.corrvec) - std(Pstat_act.LC.SC.corrvec);
                end
            end

            for b = 1:varianceSteps
                act_input_var = 4 * 2^(2*(b-1));
                act_input_mean = 50 * 0.5^(2*(b-1));
                for r = 1:params.nRerunVar
                    printCounter(r, 'maxVal', params.nRerunVar, 'stringVal', sprintf('Variance step %d',b));
                    [Pstat_act,~] = simCorr(data_unit_num, corrbeta_P, 0, params.nBinShort, act_input_mean, act_input_var);
                    ff_p(b,r) = mean(Pstat_act.HC.SC.ffvec);
                    ori_dstd_varstep(b,r) = std(Pstat_act.HC.PO.SC.corrvec) - std(Pstat_act.HC.NO.SC.corrvec);
                    ctr_dstd_varstep(b,r) = std(Pstat_act.HC.SC.corrvec) - std(Pstat_act.LC.SC.corrvec);            
                end
            end
        end
        
        if ~params.bulkSim
            save('pop_save.mat','Pstat','RGstat', ...
                'delta_std_ctr_P','delta_std_ctr_RG','delta_std_ori_P', 'delta_std_ori_RG','delta_mean_ctr_P','delta_mean_ctr_RG','delta_mean_ori_P','delta_mean_ori_RG', ...
                'mp_std','ori_dstd','ctr_dstd','ff_p','ori_dstd_varstep','ctr_dstd_varstep','simSmallBin');
        end
    else
        load('pop_save.mat')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [ecker_sc_lc, ecker_sc_hc] = separateContrast(all_spikes, conditions);    
    [cEsc_lc,mEsc_lc,~] = getStats(ecker_sc_lc);
    [cEsc_hc,mEsc_hc,~] = getStats(ecker_sc_hc);          
    
    [~,~,~,~,cEsc_po,cEsc_no] = separateOrientation(all_spikes, conditions);
    
    delta_std_ctr_E_all = std(cEsc_hc)-std(cEsc_lc);
    delta_std_ori_E_all = std(cEsc_po)-std(cEsc_no);
    delta_mean_ctr_E_all = mean(cEsc_hc)-mean(cEsc_lc);
    delta_mean_ori_E_all = mean(cEsc_po)-mean(cEsc_no);
    
    rand(10); % this is here to make the random stream step

    delta_std_ctr_E = zeros(params.nBootstrap,1);
    delta_std_ori_E = zeros(params.nBootstrap,1);
    delta_mean_ctr_E = zeros(params.nBootstrap,1);
    delta_mean_ori_E = zeros(params.nBootstrap,1);
    
    cEsc_hc_sh = cEsc_hc(randperm(length(cEsc_hc)));
    cEsc_lc_sh = cEsc_lc(randperm(length(cEsc_lc)));
    cEsc_po_sh = cEsc_po(randperm(length(cEsc_po)));
    cEsc_no_sh = cEsc_no(randperm(length(cEsc_no)));
    
    bootsize_hc = floor(length(cEsc_hc)/params.nBootstrap);
    bootsize_lc = floor(length(cEsc_lc)/params.nBootstrap);
    bootsize_po = floor(length(cEsc_po)/params.nBootstrap);
    bootsize_no = floor(length(cEsc_no)/params.nBootstrap);
    for b = 1:params.nBootstrap            
        endindex_hc = b*bootsize_hc;
        endindex_lc = b*bootsize_lc;
        endindex_po = b*bootsize_po;
        endindex_no = b*bootsize_no;
        if b==params.nBootstrap
            endindex_hc = length(cEsc_hc);
            endindex_lc = length(cEsc_lc);
            endindex_po = length(cEsc_po);
            endindex_no = length(cEsc_no);
        end
        cEsc_hc_act = cEsc_hc_sh((b-1)*bootsize_hc + 1 : endindex_hc);
        cEsc_lc_act = cEsc_lc_sh((b-1)*bootsize_lc + 1 : endindex_lc);
        cEsc_po_act = cEsc_po_sh((b-1)*bootsize_po + 1 : endindex_po);
        cEsc_no_act = cEsc_no_sh((b-1)*bootsize_no + 1 : endindex_no);
                
        delta_std_ctr_E(b) = std(cEsc_hc_act)-std(cEsc_lc_act);
        delta_std_ori_E(b) = std(cEsc_po_act)-std(cEsc_no_act);
        delta_mean_ctr_E(b) = mean(cEsc_hc_act)-mean(cEsc_lc_act);
        delta_mean_ori_E(b) = mean(cEsc_po_act)-mean(cEsc_no_act);
    end    
    
%     [delta_std_ctr_E,delta_std_ori_E,delta_mean_ctr_E,delta_mean_ori_E] = histogramStability(false,all_spikes,conditions);
    
    e_dstd_ori = std(cEsc_po) - std(cEsc_no);
    e_dstd_ctr = std(cEsc_hc) - std(cEsc_lc);    
    
    [~,pVals(1),~,~] = ttest(delta_std_ctr_P);
    [~,pVals(2),~,~] = ttest(delta_std_ctr_RG);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~ params.bulkSim

        fprintf('\n');
        fprintf('Data HC rate mean %.3f Hz\n',mean(mEsc_hc) / dataBinWidthSec);
        fprintf('Data LC rate mean %.3f Hz\n',mean(mEsc_lc) / dataBinWidthSec);
        
        simBinWidthSec = (params.nBinLong * simSmallBin) / 1000;
        
        fprintf('\n');
        fprintf('DS HC rate mean %.3f Hz\n',mean(Pstat.HC.SC.meanvec) / simBinWidthSec);
        fprintf('DS LC rate mean %.3f Hz\n',mean(Pstat.LC.SC.meanvec) / simBinWidthSec);        
        fprintf('DS HC ff mean %.3f\n',mean(Pstat.HC.SC.ffvec));
        fprintf('DS LC ff mean %.3f\n',mean(Pstat.LC.SC.ffvec));

        fprintf('\n');
        fprintf('RG HC rate mean %.3f Hz\n',mean(RGstat.HC.SC.meanvec)/ simBinWidthSec);
        fprintf('RG LC rate mean %.3f Hz\n',mean(RGstat.LC.SC.meanvec)/ simBinWidthSec);
        rg_ff = RGstat.HC.SC.ffvec;
        rg_ff(isnan(rg_ff)) = [];
        fprintf('RG HC ff mean %.3f\n',mean(rg_ff));                    
        rg_ff_lc = RGstat.LC.SC.ffvec;
        rg_ff_lc(isnan(rg_ff_lc)) = [];
        fprintf('RG LC ff mean %.3f\n',mean(rg_ff_lc));
        fprintf('\n');

    end
    
    if params.makePlots
        close all                
        res = 40;                
        secondRowPanels = 4;        
        grey_p = 0.1;
        grey_r = 0.7;
        red_exp = [0.7 0.2 0.1];
          
        figCtr = figure('units','normalized','outerposition',[0 0 1 1]);
        legends = {'HC', 'LC'};
        firstRowPanels = 4;
        plotHistogramRow(2,firstRowPanels,[1 2 3],{Pstat.HC.SC.corrvec Pstat.LC.SC.corrvec; RGstat.HC.SC.corrvec RGstat.LC.SC.corrvec; cEsc_hc cEsc_lc},{'Doubly stochastic Poisson','Rectified Gaussian','Experiment'},{},{},res,'spike count correlations')
        equalizePlots(2,firstRowPanels,1,1,1,3,true,true)        
        subplot(2,firstRowPanels,3)
        legend(legends);
        legend boxoff
        inset(2,firstRowPanels,1,@() doubleHistPlot(Pstat.HC.MP.corrvec, Pstat.LC.MP.corrvec, [-1 1], res, {}, [], false),{'membrane potential', 'correlations'},[-1 1])
        inset(2,firstRowPanels,2,@() doubleHistPlot(RGstat.HC.MP.corrvec, RGstat.LC.MP.corrvec, [-1 1], res, {}, [], false),{'membrane potential', 'correlations'},[-1 1])
        %inset(2,firstRowPanels,1,@() doubleHistPlot(Pstat.HC.MP.corrvec, [], [-1 1], res, {}, [], false),{'membrane potential', 'correlations'},[-1 1])
        %inset(2,firstRowPanels,2,@() doubleHistPlot(RGstat.HC.MP.corrvec, [], [-1 1], res, {}, [], false),{'membrane potential', 'correlations'},[-1 1])
          
        subplot(2,firstRowPanels,4)
        topDataAxis = plotRunHists(delta_std_ctr_P,delta_std_ctr_RG,'dataMinTrial',params.minDataTrials,'legend',true,'xlabel',true,'scattercol',red_exp);
        
        subplot(2,secondRowPanels,5)        
        %barPlot({delta_mean_ctr_P, delta_mean_ctr_RG, delta_mean_ctr_E_all},{'DSP','RG','Exp'},[],'','\Delta corr. mean (HC - LC)','zeroDiff','CI');
        [expts_ctr_mean,expps_ctr_mean] = barPlot({[0], [0], delta_mean_ctr_E},{'DSP','RG','Exp'},[],'','\Delta corr. mean (HC - LC)','zeroDiff','SEM',false,repmat(red_exp,3,1));                
        simps_ctr_mean = plotConfBars({delta_mean_ctr_P,delta_mean_ctr_RG},[grey_p grey_r],0.95);
        box off
        
        subplot(2,secondRowPanels,6)
        hold on
        %barPlot({delta_std_ctr_P, delta_std_ctr_RG, delta_std_ctr_E_all},{'DSP','RG','Exp'},[],'','\Delta corr. SD (HC - LC)','zeroDiff','CI',false,[],[true,true,true]);        
        [expts_ctr_sd,expps_ctr_sd] = barPlot({[0], [0], delta_std_ctr_E},{'DSP','RG','Exp'},[],'','\Delta SD_{ctr}','zeroDiff','SEM',false,repmat(red_exp,3,1));        
        simps_ctr_sd = plotConfBars({delta_std_ctr_P,delta_std_ctr_RG},[grey_p grey_r],0.95);
        equalizePlots(2,4,2,2,1,2,false,true,true)
        
        if params.simVariedMPCorr
            subplot(2,secondRowPanels,7)
            scatterErrorPlot(mp_std,ctr_dstd,e_dstd_ctr,params.nRerunVar,'MP corr. SD','\Delta SD_{ctr}',red_exp);

            subplot(2,secondRowPanels,8)
            scatterErrorPlot(ff_p,ctr_dstd_varstep,e_dstd_ctr,params.nRerunVar,'Mean SC Fano factor','\Delta SD_{ctr}',red_exp);
        end

        saveFigure(figCtr,'fig8_pop_contrast_data',false,'coords',{});                 
        
        figOri = figure('units','normalized','outerposition',[0 0 1 1]);
        legends = {'PO', 'NO'};
        firstRowPanels = 3;
        plotHistogramRow(2,3,[1 2 3],{Pstat.HC.PO.SC.corrvec Pstat.HC.NO.SC.corrvec; RGstat.HC.PO.SC.corrvec RGstat.HC.NO.SC.corrvec; cEsc_po cEsc_no},{'Doubly stochastic','Rectified Gaussian','Experiment'},{},{},res,'spike count correlations')
        equalizePlots(2,3,1,1,1,3,true,true)
        subplot(2,3,3)
        legend(legends);        
        legend boxoff
        % TODO only plot or-pref related MP corr distributions
        inset(2,3,1,@() doubleHistPlot(Pstat.HC.MP.corrvec, Pstat.LC.MP.corrvec, [-1 1], res, {}, [], false),{'membrane potential', 'correlations'},[-1 1])
        inset(2,3,2,@() doubleHistPlot(RGstat.HC.MP.corrvec, RGstat.LC.MP.corrvec, [-1 1], res, {}, [], false),{'membrane potential', 'correlations'},[-1 1])        
                  
        subplot(2,secondRowPanels,5)        
        [expts_ori_mean,expps_ori_mean] = barPlot({[0], [0], delta_mean_ori_E},{'DSP','RG','Exp'},[],'','\Delta corr. mean (PO - NO)','zeroDiff','SEM',false,repmat(red_exp,3,1));                
        simps_ori_mean = plotConfBars({delta_mean_ori_P,delta_mean_ori_RG},[grey_p grey_r],0.95);
        box off
        
        subplot(2,secondRowPanels,6)
        hold on
        [expts_ori_sd,expps_ori_sd] = barPlot({[0], [0], delta_std_ori_E},{'DSP','RG','Exp'},[],'','\Delta SD_{ori}','zeroDiff','SEM',false,repmat(red_exp,3,1));        
        simps_ori_sd = plotConfBars({delta_std_ori_P,delta_std_ori_RG},[grey_p grey_r],0.95);
        equalizePlots(2,4,2,2,1,2,false,true,true)
        
        if params.simVariedMPCorr
            subplot(2,secondRowPanels,7)
            scatterErrorPlot(mp_std,ori_dstd,e_dstd_ori,params.nRerunVar,'MP corr. SD','\Delta SD_{ori}',red_exp);

            subplot(2,secondRowPanels,8)
            scatterErrorPlot(ff_p,ori_dstd_varstep,e_dstd_ori,params.nRerunVar,'Mean SC Fano factor','\Delta SD_{ori}',red_exp);
        end

        saveFigure(figOri,'fig10_pop_orient_data',false);     
    end

%     ressub.nUnits = n_unit;
%     ressub.nRerun = params.nRerun;    
%     
%     ressub.P_MP_corr_sd = std(Pstat.HC.MP.corrvec);
%     ressub.RG_MP_corr_sd = std(RGstat.HC.MP.corrvec);
%     ressub.MP_mean_factor = params.contrast_mean_factor;
%     
%     ressub.P_HC_rate_mean = mean(Pstat.HC.SC.meanvec)/ simBinWidthSec;
%     ressub.P_HC_ff_mean = mean(Pstat.HC.SC.ffvec); 
%     ressub.P_HC_corr_mean = mean(Pstat.HC.SC.corrvec); 
%     ressub.P_LC_corr_mean = mean(Pstat.LC.SC.corrvec); 
%     ressub.P_HC_corr_sd = std(Pstat.HC.SC.corrvec); 
%     ressub.P_LC_corr_sd = std(Pstat.LC.SC.corrvec); 
%     
%     ressub.P_PO_corr_mean = mean(Pstat.HC.PO.SC.corrvec); 
%     ressub.P_NO_corr_mean = mean(Pstat.HC.NO.SC.corrvec); 
%     ressub.P_PO_corr_sd = std(Pstat.HC.PO.SC.corrvec); 
%     ressub.P_NO_corr_sd = std(Pstat.HC.NO.SC.corrvec); 
%     
%     ressub.P_ctrdsd_p = simps_ctr_sd(1);
%     ressub.P_ctrdmean_p = simps_ctr_mean(1);
%     ressub.P_oridsd_p = simps_ori_sd(1);
%     ressub.P_oridmean_p = simps_ori_mean(1);
%     
%     ressub.RG_HC_rate_mean = mean(RGstat.HC.SC.meanvec)/ simBinWidthSec;    
%     ressub.RG_HC_ff_mean = mean(RGstat.HC.SC.ffvec);
%     ressub.RG_HC_corr_mean = mean(RGstat.HC.SC.corrvec); 
%     ressub.RG_LC_corr_mean = mean(RGstat.LC.SC.corrvec); 
%     ressub.RG_HC_corr_sd = std(RGstat.HC.SC.corrvec); 
%     ressub.RG_LC_corr_sd = std(RGstat.LC.SC.corrvec); 
%     
%     ressub.RG_PO_corr_mean = mean(RGstat.HC.PO.SC.corrvec); 
%     ressub.RG_NO_corr_mean = mean(RGstat.HC.NO.SC.corrvec); 
%     ressub.RG_PO_corr_sd = std(RGstat.HC.PO.SC.corrvec); 
%     ressub.RG_NO_corr_sd = std(RGstat.HC.NO.SC.corrvec); 
%     
%     ressub.RG_ctrdsd_p = simps_ctr_sd(2);
%     ressub.RG_ctrdmean_p = simps_ctr_mean(2);
%     ressub.RG_oridsd_p = simps_ori_sd(2);
%     ressub.RG_oridmean_p = simps_ori_mean(2);
%     
%     ressub.EXP_HC_corr_mean = mean(cEsc_hc);
%     ressub.EXP_LC_corr_mean = mean(cEsc_lc);
%     ressub.EXP_HC_corr_sd = std(cEsc_hc);
%     ressub.EXP_LC_corr_sd = std(cEsc_lc);            
%     ressub.EXP_PO_corr_mean = mean(cEsc_po);
%     ressub.EXP_NO_corr_mean = mean(cEsc_no);
%     ressub.EXP_PO_corr_sd = std(cEsc_po);
%     ressub.EXP_NO_corr_sd = std(cEsc_no);            
%     
%     ressub.EXP_ctrdsd_t = expts_ctr_sd(3);
%     ressub.EXP_ctrdsd_p = expps_ctr_sd(3);
%     ressub.EXP_ctrdmean_t = expts_ctr_mean(3);
%     ressub.EXP_ctrdmean_p = expps_ctr_mean(3);
%     ressub.EXP_oridsd_t = expts_ori_sd(3);
%     ressub.EXP_oridsd_p = expps_ori_sd(3);
%     ressub.EXP_oridmean_t = expts_ori_mean(3);
%     ressub.EXP_oridmean_p = expps_ori_mean(3);
%     
%     substituteText('jneurophys/results_pop_fig8_9_10.template',ressub,'##',false,'jneurophys/results_pop_fig8_9_10.tex');
end    

function pvals = plotConfBars(data,colors,conflevel)
    pvals = zeros(length(data),1);
    for i = 1:length(data)
        [~,~,pvals(i)] = plotConfinf(data{i},conflevel,i,colors(i)*ones(1,3),false,false,'none');
    end    
    for i = 1:length(data)
        plotConfinf(data{i},conflevel,i,colors(i)*ones(1,3),false,true,'star');
    end
end

function plotHistogramRow(rownum,colnum,plotindexes,data,ylabelstr,titles,legends,resolution,xlabelstr)
    % data has to be a colnum x 2 cell array
    for i = 1:length(plotindexes)        
        subplot(rownum,colnum,plotindexes(i))
        doubleHistPlot(data{i,1}, data{i,2}, [-0.3 0.6], resolution, legends, [], false)
        if ~iscell(ylabelstr)
            ylabelstr = {ylabelstr};
        end
        if i <= length(ylabelstr)
            ylabel(ylabelstr{i})
        end
        if ~isempty(titles)
            title(titles{i})
        end
        if nargin > 8
            xlabel(xlabelstr);
        end
    end
end

function inset(rnum,cnum,spidx,call,text,limits)
    subplot(rnum,cnum,spidx)
    pos = get(gca,'Position');
    inset_width = pos(3) * 0.3 + 0.01;
    inset_height = pos(4) * 0.3 + 0.01;
    x_start = pos(1) + pos(3) - inset_width - 0.01;
    y_start = pos(2) + pos(4) - inset_height - 0.01;
    axes('Position',[x_start y_start inset_width inset_height])    
    call();
    set(gca,'XTickLabel',{},'box','on');
    xlh = xlabel(text,'FontSize',12);
    xlh.Position(2) = xlh.Position(2) - 0.01;
    xlim(limits)
end

function setrandseed(randseed)
    if strcmp(randseed,'last')
        % TODO check if exsits
        load('lastrandseed');
    elseif strcmp(randseed,'leave')
        return
    end
    s = RandStream('mt19937ar','Seed',randseed);
    RandStream.setGlobalStream(s);
    randseed = s.Seed;
    save('lastrandseed.mat','randseed');
end

function scatterErrorPlot(xdata,ydata,expval,nRerun,xlab,ylab,expcolor)
    hold on
    set(gca,'FontSize',16);
    scatter(xdata(:),ydata(:),50, 0.5 * [1 1 1],'filled')
    xdataMean = mean(xdata,2);
    errorbar(xdataMean,mean(ydata,2),std(ydata,0,2) ./ sqrt(nRerun),'k','LineWidth', 2)
    xexcess = (max(xdata(:)) - min(xdata(:))) * 0.1;    
    xl = [min(xdata(:))-xexcess max(xdata(:))+xexcess];
    xlim(xl);
    plot(xl,[expval expval],'--','Color',expcolor,'LineWidth', 2)
    text(min(xdata(:)) + (max(xdata(:)) - min(xdata(:))) * 0.6,expval*0.9,'Experiment','FontSize',12)
    plot(xl,[0 0],':k','LineWidth', 1)
    
    trueymax = max(max(ydata(:)),expval);
    trueymin = min(min(ydata(:)),expval);
    yexcess = (trueymax - trueymin) * 0.1;
    yl = [trueymin-yexcess trueymax+yexcess];
    ylim(yl);
    
    for v = 1:size(xdata,1)
        [~,pval,~,~] = ttest(ydata(v,:));
        % fprintf('%s p value: %.4f t stat: %.4f df: %d\n',labels{v},pval,stats.tstat,stats.df);
        if pval < 0.05
            starstring = '*';
            % yl = ylim();
            text(xdataMean(v), yl(2) - 0.07 * (yl(2)-yl(1)), starstring, 'FontSize', 30, 'HorizontalAlignment', 'center');                
        end
    end
    
    xlim(xl)
    xlabel(xlab,'FontSize',16);
    ylabel(ylab,'FontSize',16); 
end

