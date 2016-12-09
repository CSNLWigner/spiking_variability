function [s, r, u] = gen_spikes_nonlin_poiss(sampleNo, mu, Sigma, k, Vth, m, nonLinExponent)

% USAGE: [s, r, u] = gen_spikes_nonlin_poiss(sampleNo, mu, Sigma, k, Vth, m, nonLinExponent)
%
% The function generates spikes according to a Doubly Stochastic Poisson
% model (DSP). The process assumes that a multivariate Gaussian
% subthreshold activity is transformed by the firing rate nonlinearity to
% obtain firing rates.  The firing rate nonlinearity is a threshold
% power-law nonlinearity:
%
% r = k * (V-Vth)^m_+
%
% From the instanteous firing rate spikes are obtained through a Poisson
% process.
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

N=length(mu);

u = randnorm(sampleNo(1)*sampleNo(2),mu,[],Sigma);

fu=zeros(size(u));
fu(u<0)=-abs(u(u<0)).^nonLinExponent;
fu(u>=0)=u(u>=0).^nonLinExponent;

uRect = fu - Vth; uRect(uRect<0)=0;

r = k * uRect.^m;

s=zeros(size(mu,1),prod(sampleNo(1)));
for i=1:sampleNo(1),
    s(:,i) = sum(poissrnd(r(:,(i-1)*sampleNo(2)+1:i*sampleNo(2)),N,sampleNo(2)),2);
end

u=reshape(u,size(u,1),sampleNo(2),sampleNo(1));