function [cm,actualMean] = sampleNonzeroMeanCorrmat(dim,dispersionParam,targetMean)
    if abs(targetMean) > 1
        error('E!');
    end
    if targetMean < 0
        error('Does not work well for negative expected values.')
    end
    if dim < 30
        warning('The mean of correlation matrices of dimension less than 30 will vary substantially due to sampling noise.')
    end
    A = eye(dim) + nodiag(targetMean*2 + dispersionParam * randn(dim,dim));
    B = nearestSPD(A);
    cm = corrcov(B);
    currentMean = mean(triuvals(cm));
    cm = eye(dim) + nodiag((targetMean/currentMean) * cm);
    actualMean = mean(triuvals(cm));
end