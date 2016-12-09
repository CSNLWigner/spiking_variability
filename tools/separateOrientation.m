function [mean_po, mean_no, ff_po, ff_no, corr_po, corr_no] = separateOrientation(sc_all, cond)
    cond_hc = {};
    
    mean_po = [];
    ff_po = [];
    mean_no = [];
    ff_no = [];
    corr_po = [];
    corr_no = [];
    
    for sess = 1:length(sc_all)        
        [nUnit, nCond, nTrial] = size(sc_all{sess});
        % get HC trials
        lc_idx = 1:2:nCond; % odd indices seem to belong to lower contrasts, it's not consistently coded accross sessions, sometimes 1 and 10, sometimes 10 and 100
        hc_idx = lc_idx + 1;
        nHCCond = length(hc_idx);
        sc_hc = sc_all{sess}(:,hc_idx,:);
        cond_hc{end+1} = cond{sess}(hc_idx);
                
        unit_average_sc = mean(reshape(sc_hc,[nUnit, nHCCond * nTrial]),2);
        cond_average_sc = mean(sc_hc,3); % nUnit x nHCCond
        
        pref_idx = cond_average_sc > repmat(unit_average_sc,[1 nHCCond]);
        nonpref_idx = cond_average_sc <= repmat(unit_average_sc,[1 nHCCond]);        
        
        act_mean_po = zeros([nUnit 1]);
        act_var_po = zeros([nUnit 1]);
        act_mean_no = zeros([nUnit 1]);
        act_var_no = zeros([nUnit 1]);
        act_corr_po = [];
        act_corr_no = [];
        for u = 1:nUnit
            act_pref_sc = squeeze(sc_hc(u,pref_idx(u,:),:));
            act_nonpref_sc = squeeze(sc_hc(u,nonpref_idx(u,:),:));
            act_mean_po(u) = mean(act_pref_sc(:));
            act_mean_no(u) = mean(act_nonpref_sc(:));
            act_var_po(u) = var(act_pref_sc(:));
            act_var_no(u) = var(act_nonpref_sc(:));
            for u2 = u+1:nUnit
                both_pref_idx = pref_idx(u,:) & pref_idx(u2,:);
                both_nonpref_idx = nonpref_idx(u,:) & nonpref_idx(u2,:);
                
                if sum(both_pref_idx) > 0
                    sc_pref = sc_hc([u; u2],both_pref_idx,:); 
                    pref_corr = corr(reshape(sc_pref, [2 sum(both_pref_idx) * nTrial])');
                    act_corr_po = [act_corr_po; pref_corr(1,2)];
                end
                if sum(both_nonpref_idx) > 0
                    sc_nonpref = sc_hc([u; u2],both_nonpref_idx,:); 
                    nonpref_corr = corr(reshape(sc_nonpref, [2 sum(both_nonpref_idx) * nTrial])');
                    act_corr_no = [act_corr_no; nonpref_corr(1,2)];
                end
            end
        end
        
        mean_po = [mean_po; act_mean_po];
        mean_no = [mean_no; act_mean_no];
        ff_po = [ff_po; act_var_po ./ act_mean_po];
        ff_no = [ff_no; act_var_no ./ act_mean_no];
        corr_po = [corr_po; act_corr_po];
        corr_no = [corr_no; act_corr_no];
        
    end
end
