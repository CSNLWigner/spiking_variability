function [sc_lc, sc_hc] = separateContrast(sc_all, cond)
    sc_lc = {};
    sc_hc = {};
    for sess = 1:length(sc_all)
        [nUnit, nCond, nTrial] = size(sc_all{sess});
        lc_idx = 1:2:nCond; % odd indices seem to belong to lower contrasts, it's not consistently coded accross sessions, sometimes 1 and 10, sometimes 10 and 100
        hc_idx = lc_idx + 1;        
        sc_lc{end+1} = reshape(sc_all{sess}(:,lc_idx,:),[nUnit length(lc_idx)*nTrial]);
        sc_hc{end+1} = reshape(sc_all{sess}(:,hc_idx,:),[nUnit length(hc_idx)*nTrial]); 
    end
end
