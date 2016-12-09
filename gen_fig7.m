%% Figure 7

muRGhc=0.9;
vRGhc=2.75;
muRGlc=0.29;
vRGlc=3.35;
cRG=.13;

muPhc=1.42;
vPhc=.2;
muPlc=1.06;
vPlc=0.58;
cP=.95;

k=20;
m=1.4;
trialSampleNo=5;

trialNo=100000;

[sRGhc] = gen_spikes_nonlin_nonpoiss([trialNo trialSampleNo], ones(2,1)*muRGhc, [1 cRG; cRG 1]*vRGhc, k/50, 0, m, 1);
cc=corrcoef(sRGhc');
cRGhc=cc(1,2);
mRGschc=mean(sRGhc,2);
vRGschc=var(sRGhc,[],2);

[sRGlc] = gen_spikes_nonlin_nonpoiss([trialNo trialSampleNo], ones(2,1)*muRGlc, [1 cRG; cRG 1]*vRGlc, k/50, 0, m, 1);
cc=corrcoef(sRGlc');
cRGlc=cc(1,2);
mRGsclc=mean(sRGlc,2);
vRGsclc=var(sRGlc,[],2);

[sPhc] = gen_spikes_nonlin_poiss([trialNo trialSampleNo], ones(2,1)*muPhc, [1 cP; cP 1]*vPhc, k/50, 0, m, 1);
cc=corrcoef(sPhc');
cPhc=cc(1,2);
mPschc=mean(sPhc,2);
vPschc=var(sPhc,[],2);

[sPlc] = gen_spikes_nonlin_poiss([trialNo trialSampleNo], ones(2,1)*muPlc, [1 cP; cP 1]*vPlc, k/50, 0, m, 1);
cc=corrcoef(sPlc');
cPlc=cc(1,2);
mPsclc=mean(sPlc,2);
vPsclc=var(sPlc,[],2);

%%

figID=11;
figure(figID)
clf

lw=2;

bp=1;
bw=1.4;
sping=0.6;

subplot(1,3,1)
hold on
plot([1 2], [mRGsclc(1) mRGschc(1)],'-d','color',colorRG,'markersize',12,'linewidth',lw)
plot([1 2], [mPsclc(1) mPschc(1)],'-o','color',colorP,'markersize',12,'linewidth',lw)
set(gca,'xtick',[1 2],'xticklabel',{'LC' 'HC'},'fontsize',16,...
    'ylim',[0 5],'xlim',[0.7 2.3],'ytick',0:2.5:10,'yticklabel',{0 25 50 75 100})
ylabel('firing rate')

subplot(1,3,2)
hold on
plot([1 2], [vRGsclc(1)/mRGsclc(1) vRGschc(1)/mRGschc(1)],'-d','color',colorRG,'markersize',12,'linewidth',lw)
plot([1 2], [vPsclc(1)/mPsclc(1) vPschc(1)/mPschc(1)],'-o','color',colorP,'markersize',12,'linewidth',lw)
set(gca,'xtick',[1 2],'xticklabel',{'LC' 'HC'},'fontsize',16,...
    'ylim',[0 1.5],'xlim',[0.7 2.3],'ytick',0:.5:3)
ylabel('Fano factor')

subplot(1,3,3)
hold on
plot([1 2], [cRGlc cRGhc],'-d','color',colorRG,'markersize',12,'linewidth',lw)
plot([1 2], [cPlc cPhc],'-d','color',colorP,'markersize',12,'linewidth',lw)
set(gca,'xtick',[1 2],'xticklabel',{'LC' 'HC'},'fontsize',16,...
    'ylim',[0 0.3],'xlim',[0.7 2.3],'ytick',0:.1:1)
ylabel('spike count correlation')

set(gcf,'color','white')
scrsz = get(0,'ScreenSize');
fxDim=1000;
fyDim=280;
set(gcf,'Position',[1 scrsz(4)-fyDim fxDim fyDim])