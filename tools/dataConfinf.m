function [mu,CI] = dataConfinf(x,conflevel)
    mu = mean(x);
%     SD = std(x);
%     SEM = std(x)/sqrt(length(x));               % Standard Error
%     margin = (1 - conflevel) / 2;
%     ts = tinv([margin  1-margin],length(x)-1);      % T-Score
%     CI = mu + ts*SD
    
    % TODO conflevels other than 95
    [cdf,x_eval] = ecdf(x);
    %cdf
    
    margin = (1 - conflevel) / 2;
    low_bound = x_eval(1);
    high_bound = x_eval(end);
    for i = 2:length(x_eval)
        if cdf(i) > margin && cdf(i-1) <= margin
            low_bound = x_eval(i);
        elseif cdf(i) > (1 - margin) && cdf(i-1) <= (1 - margin)
            high_bound = x_eval(i);
        end
    end
    CI = [low_bound high_bound];
end