%% Figure 4

muRGorig=0.9;
muPorig=1.5;
cRG=.13;
cP=.95;
vRG=2.75;
vP=.2;
k=20;
m=1.4;
trialSampleNo=5;

trialNo=100000;

rg=-1.5:.1:2;

rgl=numel(rg);

mPscs=zeros(2,rgl);
mRGscs=zeros(2,rgl);
cPscs=zeros(1,rgl);
cRGscs=zeros(1,rgl);
vPscs=zeros(2,rgl);
vRGscs=zeros(2,rgl);
muRGrg=zeros(1,rgl);
muPrg=zeros(1,rgl);
for i=1:rgl,
    muRG = muRGorig -.9 +rg(i);
    muP = muPorig -.55 +.6*rg(i);
    muRGrg(i)=muRG;
    muPrg(i)=muP;
    
    [sP] = gen_spikes_nonlin_poiss([trialNo trialSampleNo], ones(2,1)*muP, [1 cP; cP 1]*vP, k/50, 0, m, 1); 
    cc=corrcoef(sP');
    cPscs(:,i)=cc(1,2);
    mPscs(:,i)=mean(sP,2);
    vPscs(:,i)=var(sP,[],2);

    [sRG] = gen_spikes_nonlin_nonpoiss([trialNo trialSampleNo], ones(2,1)*muRG, [1 cRG; cRG 1]*vRG, k/50, 0, m, 1);
    cc=corrcoef(sRG');
    cRGscs(:,i)=cc(1,2);
    mRGscs(:,i)=mean(sRG,2);
    vRGscs(:,i)=var(sRG,[],2);
end


%%
figID = 5;
figure(figID)
clf

subplot(2,3,1)
hold on;
plot([0 7],[0 7],'color','black','linewidth',1)
plot(mPscs(1,:),mRGscs(1,:),'o','color',colorL)
set(gca,'xlim',[0 7],'ylim',[0 7],'fontsize',16,'box','off','dataaspectratio',[1 1 1],...
    'xtick',0:2:6,'xticklabel',{0 20 40 60},'ytick',0:2:6,'yticklabel',{0 20 40 60})
xlabel('doubly stochastic Poisson')
ylabel('rectified Gaussian')

subplot(2,3,4)
hold on
plot(muRGrg,mPscs(1,:),'color',colorP,'linewidth',2)
plot(muRGrg,mRGscs(1,:),'color',colorRG,'linewidth',2)
plot(muRGorig, (mRGscs(1,min(find(muRGrg>muRGorig-0.0001)))+mPscs(1,min(find(muRGrg>muRGorig-0.0001))))/2,'o','markersize',12,'linewidth',2,'color',ones(1,3)*0)
set(gca,'fontsize',16,'box','off',...
    'ylim',[0 7],'ytick',0:2:6,'yticklabel',{0 20 40 60},...
    'xtick',-2:2,'xticklabel',{-2 '' 0 '' 2})
% ylabel('spike count')
ylabel('firing rate')

subplot(2,3,2)
hold on;
plot([0 1.5],[0 1.5],'color','black','linewidth',1)
plot(vPscs(1,:)./mPscs(1,:),vRGscs(1,:)./mRGscs(1,:),'o','color',colorL)
set(gca,'ylim',[0 1.5],'fontsize',16,'box','off','dataaspectratio',[1 1 1],...
    'xtick',0:.5:1.5,'xticklabel',{0 '' 1 ''},'ytick',0:.5:1.5,'yticklabel',{0 '' 1 ''})
xlabel('doubly stochastic Poisson')
% ylabel('rectified Gaussian')

subplot(2,3,5)
hold on
plot(muRGrg,vRGscs(1,:)./mRGscs(1,:),'color',colorRG,'linewidth',2)
plot(muRGrg,vPscs(1,:)./mPscs(1,:),'color',colorP,'linewidth',2)
plot(muRGorig, vRGscs(1,min(find(muRGrg>muRGorig-0.0001)))/mRGscs(1,min(find(muRGrg>muRGorig-0.0001))),'o','markersize',12,'linewidth',2,'color',ones(1,3)*0)
set(gca,'ylim',[0 1.5],'fontsize',16,'box','off','ytick',0:.5:1.5,'yticklabel',{0 '' 1 ''},...
        'xtick',-2:2,'xticklabel',{-2 '' 0 '' 2})
ylabel('Fano factor')
% legend('DsP', 'RG')
% legend('boxoff')

subplot(2,3,3)
hold on;
% plot([0 .2],[0 .2],'color','black','linewidth',1)
plot([-.2 .2],[-.2 .2],'color','black','linewidth',1)
plot(cPscs(1,:),cRGscs(1,:),'o','color',colorL)
set(gca,'xlim',[0 .2],'ylim',[0 .2],'ytick',0:.1:.2,'xtick',0:.1:.2,...
    'fontsize',16,'box','off','dataaspectratio',[1 1 1])
% set(gca,'xlim',[-.2 0],'ylim',[-.2 0],'fontsize',16,'box','off','dataaspectratio',[1 1 1])
xlabel('doubly stochastic Poisson')
% ylabel('rectified Gaussian')

subplot(2,3,6)
hold on
plot(muRGrg,cPscs,'color',colorP,'linewidth',2)
plot(muRGrg,cRGscs,'color',colorRG,'linewidth',2)
plot(muRGorig, cPscs(1,min(find(muRGrg>muRGorig-0.0001))),'o','markersize',12,'linewidth',2,'color',ones(1,3)*0)
% set(gca,'ylim',[-.2 0],'fontsize',16,'box','off')
set(gca,'ylim',[0 .2],'ytick',0:.1:.2,'fontsize',16,'box','off',...
        'xtick',-2:2,'xticklabel',{-2 '' 0 '' 2})
ylabel('correlation')

set(gcf,'color','white')
scrsz = get(0,'ScreenSize');
fxDim=1000;
fyDim=650;
set(gcf,'Position',[1 scrsz(4)-fyDim fxDim fyDim])

