function [s, r, u] = gen_spikes_nonlin_nonpoiss(sampleNo, mu, Sigma, k, Vth, m, nonLinExponent)

% USAGE: [s, r, u] = gen_spikes_nonlin_nonpoiss(sampleNo, mu, Sigma, k, Vth, m, nonLinExponent)
%
% This is an implementation of the Rectified Gaussian model, where
% subthreshold activity is a result of a multivariate stochastic process
% but spikes are generated as determinsitacally as possible: the firing 
% rate determines the number of spikes up to a random initial condition for
% accumulation factor. Subthreshold activity is transformed by the firing
% rate nonlinearity to obtain firing rates.  The firing rate nonlinearity
% is a threshold power-law nonlinearity:
%
% r = k * (V-Vth)^m_+
%
% Firing rate is then integrated in order to obtain spikes: spikes are
% generated when the integral reaches integer values.
%
% INPUT PARAMETERS
% sampleNo: 1x2 vector: the first containing the number of membrane
%   potential samples to be taken, the second the length of a trial (in terms
%   of membrane potential samples in a single trial)
% mu: vector defining the mean of subthreshold activity of a
%   nuron/population of neurons
% Sigma: covariance matrix of Gaussian that defines the variance and
%   correlation structure of subthreshold activity 
% k: gain of the firing rate nonlinearity
% Vth: threshold of firing rate nonlinearity
% m: exponent of firing rate nonlinearity
% nonLineExponent: exponent of a potential nonlinearity that can be applied
% to the Gaussian distribution in order to obtain membrane potential
% distribution with a kurtosis that differs from that of the Gaussian
%
% Created by Gergo Orban


u = randnorm(sampleNo(1)*sampleNo(2),mu,[],Sigma);

fu=zeros(size(u));
fu(u<0)=-abs(u(u<0)).^nonLinExponent;
fu(u>=0)=u(u>=0).^nonLinExponent;

uRect = fu - Vth; uRect(uRect<0)=0;

r = k * uRect.^m;

s=zeros(size(mu,1),prod(sampleNo(1)));
for i=1:sampleNo(1),
    rs = sum(r(:,(i-1)*sampleNo(2)+1:i*sampleNo(2)),2);
    rr = rand(size(r,1),1);
    s(:,i) = floor(rs) + ((rs-floor(rs))>rr);
end

u=reshape(u,size(u,1),sampleNo(2),sampleNo(1));