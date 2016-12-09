function [corr,m,f] = getStats(sc)
    if ~iscell(sc)
        sc = {sc};
    end
    corr = [];
    m = [];
    f = [];
    for s=1:length(sc)
        
        summedSpikes = sum(sc{s},2); % nUnit x 1
        sc_nonzero = sc{s};
        sc_nonzero(summedSpikes == 0,:) = [];
        % cc=corrcoef(sc_nonzero');
        cc=corrcoef(zscore((sc_nonzero')));
        
        if any(isnan(cc))
            sum(sc{s},2)
            error('nan in corr');
        end
        
        corr = [corr; triuvals(cc)];
        am = mean(sc{s},2);
        m = [m; am];
        v = var(sc{s},[],2);
        f = [f; v ./ am];
    end    
end