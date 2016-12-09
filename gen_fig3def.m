%% Figure 3D-F

muRG=0.9;
muP=1.42;
% cRG=.15;
cRG=.13;
cP=.95;
vRG=2.75;
vP=.2;
k=20;
m=1.4;
trialSampleNo=5;

trialNo=200000;

[sP] = gen_spikes_nonlin_poiss([trialNo trialSampleNo], ones(2,1)*muP, [1 cP; cP 1]*vP, k/50, 0, m, 1); 
cc=corrcoef(sP');
cPsc=cc(1,2)
mPsc=mean(sP,2)
vPsc=var(sP,[],2)

[sRG] = gen_spikes_nonlin_nonpoiss([trialNo trialSampleNo], ones(2,1)*muRG, [1 cRG; cRG 1]*vRG, k/50, 0, m, 1);
cc=corrcoef(sRG');
cRGsc=cc(1,2)
mRGsc=mean(sRG,2)
vRGsc=var(sRG,[],2)

%% generation of the figure

figID=2;
figure(figID)
clf;
sNo=50;

rg=[-.5 8.5];

subplot(1,3,1)
hold on
ccP = plot_cov_contour(cov(sP'));
plot(sP(1,1:sNo)+randn(1,sNo)*0.1,sP(2,1:sNo)+randn(1,sNo)*0.1,'o','markersize',9,'linewidth',0.5,'color',ones(1,3)*0,'markerfacecolor',colorP)
plot(mean(sP(1,:)),mean(sP(2,:)),'+','markersize',12,'linewidth',1,'color',colorL,'linewidth',2)
plot(mean(sP(1,:))+ccP(:,1),mean(sP(2,:))+ccP(:,2),'linewidth',3,'color',colorL)
set(gca,'xlim',rg,'ylim', rg,'dataaspectratio',[1 1 1],'box', 'off',...
    'fontsize', 16,'xtick',[0:4:50],'ytick',[0:4:50])
% title('Doubly stochastic Poisson')
xlabel('spike count, neuron #1')
ylabel('spike count, neuron #2')

subplot(1,3,2)
ccRG = plot_cov_contour(cov(sRG'));
hold on
plot(sRG(1,1:sNo)+randn(1,sNo)*0.1,sRG(2,1:sNo)+randn(1,sNo)*0.1,'o','markersize',9,'linewidth',0.5,'color',ones(1,3)*0,'markerfacecolor',colorRG)
plot(mean(sRG(1,:)),mean(sRG(2,:)),'+','markersize',12,'linewidth',1,'color',colorL,'linewidth',2)
plot(mean(sRG(1,:))+ccRG(:,1),mean(sRG(2,:))+ccRG(:,2),'linewidth',3,'color',colorL)
set(gca,'xlim',rg,'ylim', rg,'dataaspectratio',[1 1 1],'box', 'off',...
    'fontsize', 16,'xtick',[0:4:50],'ytick',[0:4:50])
% title('Gaussian rectifier')
xlabel('spike count, neuron #1')

% set(gcf,'color','white')
% scrsz = get(0,'ScreenSize');
% fxDim=800;
% fyDim=400;
% set(gcf,'Position',[1 scrsz(4)-fyDim fxDim fyDim])

bp=1;
bw=1.4;
sping=0.6;

subplot(1,27,[20 21])
hold on
rectangle('Position',[bp-bw/2 0 bw mPsc(1)],'facecolor',colorP)
rectangle('Position',[bp-bw/2 0 bw mRGsc(1)]+[1 0 0 0]*(bw+sping),'facecolor',colorRG)
set(gca,'xlim',[0 4],'xtick',{},'fontsize',16,'ylim',[0 6],'ytick',0:2:10,'yticklabel',{0 20 40 60 80 100})
ylabel('firing rate')

subplot(1,27,[23 24])
hold on
% rectangle('Position',[bp-bw/2 0 bw vPsc(1)],'facecolor',colorP)
rectangle('Position',[bp-bw/2 0 bw vPsc(1)/mPsc(1)],'facecolor',colorP)
rectangle('Position',[bp-bw/2 0 bw vRGsc(1)/mRGsc(1)]+[1 0 0 0]*(bw+sping),'facecolor',colorRG)
set(gca,'xlim',[0 4],'xtick',{},'fontsize',16,'ytick',0:.5:1.5,'ylim',[0 1.5])
% ylabel('spike variance')
ylabel('Fano factor')

subplot(1,27,[26 27])
hold on
rectangle('Position',[bp-bw/2 0 bw cPsc(1)],'facecolor',colorP)
rectangle('Position',[bp-bw/2 0 bw cRGsc(1)]+[1 0 0 0]*(bw+sping),'facecolor',colorRG)
set(gca,'xlim',[0 4],'xtick',{},'fontsize',16,'ytick',0:.1:.5,'ylim',[0 .3])
ylabel('spike count correlation')

set(gcf,'color','white')
scrsz = get(0,'ScreenSize');
fxDim=1500;
fyDim=400;
set(gcf,'Position',[1 scrsz(4)-fyDim fxDim fyDim])
