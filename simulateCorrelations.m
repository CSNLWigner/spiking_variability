function [Pstat,RGstat,bin_length] = simulateCorrelations(n_unit, corrparam_P, corrparam_RG, n_bin, varargin)

    parser = inputParser;
    % addParameter(parser,'tc_coeff_P',50,@isnumeric); % ORIG
    % addParameter(parser,'tc_coeff_RG',100,@isnumeric); % ORIG
    
    addParameter(parser,'mp_var_b_P',4,@isnumeric);
    addParameter(parser,'mp_var_b_RG',4,@isnumeric);    
    addParameter(parser,'tc_coeff_P',50,@isnumeric);    
    addParameter(parser,'tc_coeff_RG',80,@isnumeric);
    
    addParameter(parser,'contrast_corr_factor',1,@isnumeric);    
    addParameter(parser,'contrast_mean_factor',2,@isnumeric);    
    addParameter(parser,'contrast_var_factor',0.7,@isnumeric);    
    
    addParameter(parser,'mp_exponent',1,@isnumeric);
    addParameter(parser,'rate_exponent',1.4,@isnumeric);
    
    addParameter(parser,'n_trial',1000,@isnumeric);
    
    parse(parser,varargin{:});
    params = parser.Results;

    % base_rate = 13; % Hz % ORIG
    base_rate = 7; % Hz
    theta = 0;       
    bin_length = 20; % ms    
    k = base_rate / (1000 / bin_length); % base number of spikes within a bin
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    contrast_mean_factor_P = params.contrast_mean_factor;     % from lower to higher
    contrast_var_factor_P = params.contrast_var_factor;    % from lower to higher
    
    contrast_mean_factor_RG = params.contrast_mean_factor;
    % contrast_var_factor_RG = 0.5; % ORIG
    contrast_var_factor_RG = params.contrast_var_factor;
    
    orient_range = 180;
    preferred_orientations = orient_range * rand([n_unit 1]);
    tuning_width = orient_range * 0.2;    
    preference_threshold = orient_range * 0.5;
    
    % stimulus_orientation = orient_range/2 * rand()
    stimulus_orientation = 5;
    
    pref_units = abs(preferred_orientations - stimulus_orientation) < preference_threshold;
    nonpref_units = abs(preferred_orientations - stimulus_orientation) >= preference_threshold; 
 
    mp_mean_mean_P = params.tc_coeff_P * max( normpdf(ones(n_unit,1) * stimulus_orientation, preferred_orientations, ones(n_unit,1) * tuning_width), ...
        normpdf(ones(n_unit,1) * (stimulus_orientation + orient_range/2) , preferred_orientations, ones(n_unit,1) * tuning_width));
    mp_mean_mean_RG = params.tc_coeff_RG * max( normpdf(ones(n_unit,1) * stimulus_orientation, preferred_orientations, ones(n_unit,1) * tuning_width), ...
        normpdf(ones(n_unit,1) * (stimulus_orientation + orient_range/2) , preferred_orientations, ones(n_unit,1) * tuning_width));
    
    mp_mean_std_P = 0.01;
    mp_mean_std_RG = 0.01;
    
    %mp_var_a_P = 3;    % ORIG
    mp_var_a_P = 4.5;
    %mp_var_a_RG = 3;    % ORIG
    mp_var_a_RG = 1.9;    
    %mp_var_a_RG = 1.4;    
    
    correlation_shift_dispersion = 0.05;
            
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    if length(corrparam_P(:)) == 1 
        corrmat_P = vineBeta(n_unit, corrparam_P);
    else
        % TODO test dimensionality
        corrmat_P = corrparam_P;
    end
    
%     if params.contrast_corr_factor <= 1
%         corrmat_P_LC = corrmat_P;        
%         corrmat_P_HC = offDiagMltpl(corrmat_P,params.contrast_corr_factor);
%     else
%         corrmat_P_HC = corrmat_P;
%         corrmat_P_LC = offDiagMltpl(corrmat_P,1 / params.contrast_corr_factor);
%     end        
        
    if params.contrast_corr_factor < 1
        corrmat_P_LC = corrmat_P;
        corrmat_P_HC = corrmat_P.* sampleNonzeroMeanCorrmat(n_unit,correlation_shift_dispersion,params.contrast_corr_factor);
    elseif params.contrast_corr_factor > 1        
        corrmat_P_LC = corrmat_P .* sampleNonzeroMeanCorrmat(n_unit,correlation_shift_dispersion,1 / params.contrast_corr_factor);
        corrmat_P_HC = corrmat_P;
    else
        corrmat_P_LC = corrmat_P;        
        corrmat_P_HC = corrmat_P;        
    end                
    
    mu_lc_P = mp_mean_mean_P + mp_mean_std_P * randn([n_unit 1]);
    var_lc_P = inverse_gamma(mp_var_a_P, params.mp_var_b_P, n_unit);
    C_lc_P = diag(sqrt(var_lc_P)) * corrmat_P_LC * diag(sqrt(var_lc_P));    
    
    mu_hc_P = mp_mean_mean_P * contrast_mean_factor_P + mp_mean_std_P * randn([n_unit 1]);
    var_hc_P = inverse_gamma(mp_var_a_P, params.mp_var_b_P * contrast_var_factor_P, n_unit);
    C_hc_P = diag(sqrt(var_hc_P)) * corrmat_P_HC * diag(sqrt(var_hc_P));                    
    
    Pstat.HC.MP.corrmat = corrmat_P_HC;
    Pstat.HC.MP.corrvec = triuvals(corrmat_P_HC);
    Pstat.HC.MP.meanvec = mu_hc_P;
    Pstat.HC.MP.varvec = var_hc_P;    
    Pstat.LC.MP.corrmat = corrmat_P_LC;
    Pstat.LC.MP.corrvec = triuvals(corrmat_P_LC);
    Pstat.LC.MP.meanvec = mu_lc_P;
    Pstat.LC.MP.varvec = var_lc_P;
    
    rates_P_lc = gen_rates(params.n_trial, n_bin, mu_lc_P, C_lc_P, params.mp_exponent, theta, params.rate_exponent, k);
    rates_P_hc = gen_rates(params.n_trial, n_bin, mu_hc_P, C_hc_P, params.mp_exponent, theta, params.rate_exponent, k);
    samples_P_lc = spikes_poisson(rates_P_lc);
    samples_P_hc = spikes_poisson(rates_P_hc);
    [Pstat.LC.SC.corrvec, Pstat.LC.SC.meanvec, Pstat.LC.SC.ffvec] = getStats(samples_P_lc);
    [Pstat.HC.SC.corrvec, Pstat.HC.SC.meanvec, Pstat.HC.SC.ffvec] = getStats(samples_P_hc);
    spikes_hc_po_P = samples_P_hc(pref_units, :);
    spikes_hc_no_P = samples_P_hc(nonpref_units, :);
    [Pstat.HC.PO.SC.corrvec, Pstat.HC.PO.SC.meanvec, Pstat.HC.PO.SC.ffvec] = getStats(spikes_hc_po_P);
    [Pstat.HC.NO.SC.corrvec, Pstat.HC.NO.SC.meanvec, Pstat.HC.NO.SC.ffvec] = getStats(spikes_hc_no_P);        
    
    if corrparam_RG ~= 0
        if length(corrparam_RG(:)) == 1 
            corrmat_RG = vineBeta(n_unit, corrparam_RG);
        else
            % TODO test dimensionality
            corrmat_RG = corrparam_RG;
        end
%         %corrmat_RG = vineBeta(n_unit, corrbeta_RG);
%         if params.contrast_corr_factor_RG <= 1
%             corrmat_RG_LC = corrmat_RG;
%             corrmat_RG_HC = offDiagMltpl(corrmat_RG,params.contrast_corr_factor_RG);
%         else
%             corrmat_RG_HC = corrmat_RG;
%             corrmat_RG_LC = offDiagMltpl(corrmat_RG,1 / params.contrast_corr_factor_RG);
%         end
        
        if params.contrast_corr_factor < 1
            corrmat_RG_LC = corrmat_RG;
            corrmat_RG_HC = corrmat_RG.* sampleNonzeroMeanCorrmat(n_unit,correlation_shift_dispersion,params.contrast_corr_factor);
        elseif params.contrast_corr_factor > 1        
            corrmat_RG_LC = corrmat_RG .* sampleNonzeroMeanCorrmat(n_unit,correlation_shift_dispersion,1 / params.contrast_corr_factor);
            corrmat_RG_HC = corrmat_RG;
        else
            corrmat_RG_LC = corrmat_RG;        
            corrmat_RG_HC = corrmat_RG;        
        end                      
        
        mu_lc_RG = mp_mean_mean_RG + mp_mean_std_RG * randn([n_unit 1]);
        var_lc_RG = inverse_gamma(mp_var_a_RG, params.mp_var_b_RG, n_unit);
        C_lc_RG = diag(sqrt(var_lc_RG)) * corrmat_RG_LC * diag(sqrt(var_lc_RG));

        mu_hc_RG = mp_mean_mean_RG * contrast_mean_factor_RG + mp_mean_std_RG * randn([n_unit 1]);
        var_hc_RG = inverse_gamma(mp_var_a_RG, params.mp_var_b_RG * contrast_var_factor_RG, n_unit);
        C_hc_RG = diag(sqrt(var_hc_RG)) * corrmat_RG_HC * diag(sqrt(var_hc_RG));
    
        RGstat.HC.MP.corrmat = corrmat_RG_HC;
        RGstat.HC.MP.corrvec = triuvals(corrmat_RG_HC);
        RGstat.HC.MP.meanvec = mu_hc_RG;
        RGstat.HC.MP.varvec = var_hc_RG;    
        RGstat.LC.MP.corrmat = corrmat_RG_LC;
        RGstat.LC.MP.corrvec = triuvals(corrmat_RG_LC);
        RGstat.LC.MP.meanvec = mu_lc_RG;
        RGstat.LC.MP.varvec = var_lc_RG;
        
        rates_RG_lc = gen_rates(params.n_trial, n_bin, mu_lc_RG, C_lc_RG, params.mp_exponent, theta, params.rate_exponent, k);
        rates_RG_hc = gen_rates(params.n_trial, n_bin, mu_hc_RG, C_hc_RG, params.mp_exponent, theta, params.rate_exponent, k);    
        samples_RG_lc = spikes_determ(rates_RG_lc);    
        samples_RG_hc = spikes_determ(rates_RG_hc);        
        [RGstat.LC.SC.corrvec, RGstat.LC.SC.meanvec, RGstat.LC.SC.ffvec] = getStats(samples_RG_lc);
        [RGstat.HC.SC.corrvec, RGstat.HC.SC.meanvec, RGstat.HC.SC.ffvec] = getStats(samples_RG_hc);                   
        spikes_hc_po_RG = samples_RG_hc(pref_units, :);
        spikes_hc_no_RG = samples_RG_hc(nonpref_units, :);        
        [RGstat.HC.PO.SC.corrvec, RGstat.HC.PO.SC.meanvec, RGstat.HC.PO.SC.ffvec] = getStats(spikes_hc_po_RG);
        [RGstat.HC.NO.SC.corrvec, RGstat.HC.NO.SC.meanvec, RGstat.HC.NO.SC.ffvec] = getStats(spikes_hc_no_RG);      
    else      
        RGstat = [];
    end                
end

function y = inverse_gamma(a, b, n)
    x=gamrnd(a, 1/b, [n 1]);
    y=1./x;
end

function r = gen_rates(n_trial, n_bin, mu, C, mp_exponent, Vth, rate_exponent, k)
    n_unit = length(mu);

%     u = zeros(n_unit, n_trial, n_bin);
%     for t = 1:n_trial
%         u(:,t,:) = mvnrnd(repmat(mu',[n_bin 1]), C)';
%     end
    
    u = reshape(mvnrnd(repmat(mu',[n_bin*n_trial 1]), C)',[n_unit,n_trial,n_bin]);
    
    u = (abs(u) .^ mp_exponent) .* sign(u);
    uRect = u - Vth; 
    uRect(uRect<0) = 0;
    r = k * uRect.^rate_exponent;
end

function s = spikes_poisson(r)
    s = poissrnd(r); % n_unit x n_trial x n_bin
    s = squeeze(sum(s,3)); % n_unit x n_trial
end

function s = spikes_determ(r)
    s = floor(r);
    s = squeeze(sum(s,3)); % n_unit x n_trial
end

function M2 = offDiagMltpl(M,s)
    M2 = M .* ((ones(size(M)) - eye(size(M))) * s + eye(size(M)));
    %M2 = max(min(M2,ones(size(M2))),-1*ones(size(M2)));
end