function [sc_all, conditions, trialNums, contrastLevels, unitNums, binWidth] = loadEckerData(minTrial,varargin)
    parser = inputParser;
    addParameter(parser,'evokedStartBin',30,@isnumeric);
    addParameter(parser,'evokedEndBin',70,@isnumeric);
    addParameter(parser,'contaminationThreshold',0.05,@isnumeric);
    addParameter(parser,'rateThresholdDivider',2000,@isnumeric);
    parse(parser,varargin{:});
    params = parser.Results;

    load('ecker_data_v1_binned_static.mat')
    
    sc_all = {};
    sc_all_binned = {};
    conditions = {};
    trialNums = [];
    unitNums = [];
    contrastLevels = [];
    
    binWidth = ((params.evokedEndBin - params.evokedStartBin) * 10) / 1000; % data is binned originally at 10 ms
    
    nSess = length(data);
    for sess = 1:nSess
        if size(data{sess}.spikes,4) < minTrial
            continue
        end
        
        [nUnit, nCond, nBin, nTrial] = size(data{sess}.spikes);

        % throwing out contaminated clusters
        nUnitNonCont = sum(data{sess}.contamination < params.contaminationThreshold);
        sptrEvokedNonCont = data{sess}.spikes(data{sess}.contamination < params.contaminationThreshold,:,params.evokedStartBin:params.evokedEndBin,:);

        % throwing out cells with very low rates in either condition
        minCount = ceil((params.evokedEndBin - params.evokedStartBin) * nTrial / params.rateThresholdDivider);
        badRate = [];
        for c = 1:nCond
            actEvokedSpikeTrains = sptrEvokedNonCont(:,c,:,:); % nUnitNC x 1 x nEvokedBin x nTrial
            actEvokedSpikeCount = squeeze(sum(actEvokedSpikeTrains,3)); % nUnit x nTrial
            badRate = union(badRate, find(sum(actEvokedSpikeCount,2) < minCount));
        end
        goodRateIdx = setdiff(1:nUnitNonCont,badRate);
        sptrEvokedGood = sptrEvokedNonCont(goodRateIdx,:,:,:);
        nUnitGood = length(goodRateIdx);

        sc_all{end+1} = squeeze(sum(sptrEvokedGood,3));  % nUnitGood x nCond x nTrial
        conditions{end+1} = data{sess}.conditions;
        sc_all_binned{end+1} = sptrEvokedGood;
        trialNums = [trialNums; nTrial];
        unitNums = [unitNums; nUnitGood];
        contrastLevels = [contrastLevels; data{sess}.conditions(1).contrast data{sess}.conditions(2).contrast];
    end
end
