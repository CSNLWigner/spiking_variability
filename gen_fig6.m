%% Figure 6

muRGorig=0.9;
muPorig=1.42;
muRGhc=1.2;
muRGlc=0.5;
cRG=.13;
cP=.95;
vRG=2.75;
vP=.2;
k=20;
m=1.4;
trialSampleNo=5;

trialNo=100000;

rg=[muRGhc muRGlc];
rg2=vRG:.2:vRG+1.6;
rg22=vP:.15:vP+1.2;

rgl=numel(rg);
rgl2=numel(rg2);

mPscs=zeros(2,rgl);
mRGscs=zeros(2,rgl);
cPscs=zeros(1,rgl);
cRGscs=zeros(1,rgl);
vPscs=zeros(2,rgl);
vRGscs=zeros(2,rgl);
muRGrg=zeros(1,rgl);
muPrg=zeros(1,rgl);
for i=1:rgl,
    muRG = rg(i);
    muP = .8 + .7*rg(i);
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

% muRG = muRGrg(2);
% muP = muPrg(2);
muRG = muRGorig;
muP = muPorig;

mPscs2=zeros(2,rgl2);
mRGscs2=zeros(2,rgl2);
cPscs2=zeros(1,rgl2);
cRGscs2=zeros(1,rgl2);
vPscs2=zeros(2,rgl2);
vRGscs2=zeros(2,rgl2);
vRGrg=zeros(1,rgl2);
vPrg=zeros(1,rgl2);
for i=1:rgl2,
    vRG = rg2(i);
    vP = rg22(i);
    vRGrg(i)=vRG;
    vPrg(i)=vP;
    
    [sP] = gen_spikes_nonlin_poiss([trialNo trialSampleNo], ones(2,1)*muP, [1 cP; cP 1]*vP, k/50, 0, m, 1); 
    cc=corrcoef(sP');
    cPscs2(:,i)=cc(1,2);
    mPscs2(:,i)=mean(sP,2);
    vPscs2(:,i)=var(sP,[],2);

    [sRG] = gen_spikes_nonlin_nonpoiss([trialNo trialSampleNo], ones(2,1)*muRG, [1 cRG; cRG 1]*vRG, k/50, 0, m, 1);
    cc=corrcoef(sRG');
    cRGscs2(:,i)=cc(1,2);
    mRGscs2(:,i)=mean(sRG,2);
    vRGscs2(:,i)=var(sRG,[],2);
end

%% generating the figure

figID = 9;
figure(figID)
clf

subplot(2,3,1)
hold on;
plot([0 7],[0 7],'color','black','linewidth',1)
plot(mPscs2(1,:),mRGscs2(1,:),'o','markersize',9,'color',colorL)
set(gca,'xlim',[0 7],'ylim',[0 7],'fontsize',16,'box','off','dataaspectratio',[1 1 1],...
    'xtick',[0 2 4 6],'xticklabel',{0 20 40 60},'ytick',[0 2 4 6],'yticklabel',{0 20 40 60})
xlabel('doubly stochastic Poisson')
ylabel('rectified Gaussian')

subplot(2,3,4)
hold on
plot(vRGrg,mPscs2(1,:),'-','linewidth',2,'color',colorP)
plot(vRGrg,mRGscs2(1,:),'-','linewidth',2,'color',colorRG)
plot(vRGrg(1),mRGscs2(1,1),'o','markersize',12,'color',ones(1,3)*0,'linewidth',2)
set(gca,'fontsize',16,'box','off',...
    'ylim',[0 7],'ytick',[0 2 4 6],'yticklabel',{0 20 40 60},...
    'xlim',[2.7 4.4])
ylabel('firing rate')
% legend('DS', 'RG')
% legend('boxoff')

subplot(2,3,2)
hold on;
plot([0 2],[0 2],'color','black','linewidth',1)
plot(vPscs2(1,:)./mPscs2(1,:),vRGscs2(1,:)./mRGscs2(1,:),'o','markersize',9,'color',colorL)
set(gca,'ylim',[0 1.7],'xlim',[0 1.7],'fontsize',16,'box','off',...
    'xtick',0:.5:1.5,'ytick',0:.5:1.5,'dataaspectratio',[1 1 1])
xlabel('doubly stochastic Poisson')
% ylabel('rectified Gaussian')

subplot(2,3,5)
hold on
plot(vRGrg,vPscs2(1,:)./mPscs2(1,:),'-','linewidth',2,'color',colorP)
plot(vRGrg,vRGscs2(1,:)./mRGscs2(1,:),'-','linewidth',2,'color',colorRG)
plot(vRGrg(1),vRGscs2(1,1)/mRGscs2(1,1),'o','markersize',12,'color',ones(1,3)*0,'linewidth',2)
set(gca,'fontsize',16,'box','off',...
    'ylim',[0 1.7],'ytick',[0:.5:5],...
    'xlim',[2.7 4.4])
ylabel('Fano factor')

subplot(2,3,3)
hold on;
plot([-1 1],[-1 1],'color','black','linewidth',1)
plot(cPscs2(1,:),cRGscs2(1,:),'o','markersize',9,'color',colorL)
set(gca,'xlim',[0 .41],'ylim',[0 .41],'fontsize',16,'box','off',...
    'xtick',0:.1:1,'ytick',0:.1:1,'dataaspectratio',[1 1 1])
% set(gca,'xlim',[-.2 0],'ylim',[-.2 0],'fontsize',16,'box','off','dataaspectratio',[1 1 1])
xlabel('doubly stochastic Poisson')
% ylabel('rectified Gaussian')

subplot(2,3,6)
hold on
plot(vRGrg,cPscs2,'-','linewidth',2,'color',colorP)
plot(vRGrg,cRGscs2,'-','linewidth',2,'color',colorRG)
plot(vRGrg(1),cRGscs2(1,1),'o','markersize',12,'color',ones(1,3)*0,'linewidth',2)
% set(gca,'ylim',[-.2 0],'fontsize',16,'box','off')
set(gca,'fontsize',16,'box','off',...
    'ylim',[0 .41],'ytick',0:.1:1,...
    'xlim',[2.7 4.4])
ylabel('spike count correlation')

set(gcf,'color','white')
scrsz = get(0,'ScreenSize');
fxDim=1000;
fyDim=650;
set(gcf,'Position',[1 scrsz(4)-fyDim fxDim fyDim])
