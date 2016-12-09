%% Figure 3A-C

muRG=0.9;
% muP=1.2;
muP=1.4;
% cRG=.15;
cRG=.13;
% cP=.9;
cP=.95;
vP=.2;
vRG=2.3;
k=20;
% k=10;
m=1.4;
% m=1.3;
trialSampleNo=5;

vRG_Rg=2:0.1:3;
% vRG_Rg=3:0.1:4;
muP_Rg=0.8:0.1:1.7;
cP_Rg=0.15:0.05:0.9;
vP_Rg=0.025:0.025:0.3;

rg1=vRG_Rg;
rg2=muP_Rg;
% rg1=cP_Rg;
% rg1=vP_Rg;

trialNo=30000;
% trialNo=10000;

rgl1=numel(rg1);
rgl2=numel(rg2);

mPscs=zeros(2,rgl1,rgl2);
mRGscs=zeros(2,rgl1,rgl2);
cPscs=zeros(1,rgl1,rgl2);
cRGscs=zeros(1,rgl1,rgl2);
vPscs=zeros(2,rgl1,rgl2);
vRGscs=zeros(2,rgl1,rgl2);
for i=1:rgl1,
    for j=1:rgl2,
        vRG=rg1(i);
        muP=rg2(j);
    %     cP=rg1(i);
    %     vP=rg1(i);

        [sP] = gen_spikes_nonlin_poiss([trialNo trialSampleNo], ones(2,1)*muP, [1 cP; cP 1]*vP, k/50, 0, m, 1); 
        cc=corrcoef(sP');
        cPscs(:,i,j)=cc(1,2);
        mPscs(:,i,j)=mean(sP,2);
        vPscs(:,i,j)=var(sP,[],2);

        [sRG] = gen_spikes_nonlin_nonpoiss([trialNo trialSampleNo], ones(2,1)*muRG, [1 cRG; cRG 1]*vRG, k/50, 0, m, 1);
        cc=corrcoef(sRG');
        cRGscs(:,i,j)=cc(1,2);
        mRGscs(:,i,j)=mean(sRG,2);
        vRGscs(:,i,j)=var(sRG,[],2);
    end
end

%%
figID=4;
figure(figID)
clf

subplot(1,3,1)
hold on
% plot(rg1,mPscs(1,:),rg1,mRGscs(1,:))
surf(rg2,rg1,squeeze(mPscs(1,:,:)),'facecolor',colorP)
surf(rg2,rg1,squeeze(mRGscs(1,:,:)),'facecolor',colorRG)
% zlabel('spike count')
zlabel('firing rate')
set(gca,'zlim',[0 5],'view',[60 50],'fontsize', 16,...
    'xlim',[min(rg2) max(rg2)], 'ylim',[min(rg1) max(rg1)],...
    'ztick',0:2.5:5,'zticklabel',{0 25 50},'xtick',.8:.2:1.6,'xticklabel',{'' 1, '' 1.4 ''})
% xlabel('MP mean, P')
% ylabel('spike variance, RG')

subplot(1,3,2)
hold on
surf(rg2,rg1,squeeze(vPscs(1,:,:)./mPscs(1,:,:)),'facecolor',colorP)
surf(rg2,rg1,squeeze(vRGscs(1,:,:)./mRGscs(1,:,:)),'facecolor',colorRG)
zlabel('Fano factor')
set(gca,'zlim',[0.6 1.4],'view',[60 50],'fontsize', 16,...
    'xlim',[min(rg2) max(rg2)], 'ylim',[min(rg1) max(rg1)],'ztick',.6:.4:1.4,...
    'xtick',.8:.2:1.6,'xticklabel',{'' 1, '' 1.4 ''})
xlabel('\mu_{MP(DsP)}')
ylabel('\sigma^2_{MP(RG)}')
% title(sprintf('response statistics @m=%1.1f',m))

subplot(1,3,3)
hold on
surf(rg2,rg1,squeeze(cPscs(1,:,:)),'facecolor',colorP)
surf(rg2,rg1,squeeze(cRGscs(1,:,:)),'facecolor',colorRG)
zlabel('correlation')
set(gca,'zlim',[0 0.2],'view',[60 50],'fontsize', 16,...
    'xlim',[min(rg2) max(rg2)], 'ylim',[min(rg1) max(rg1)],'ztick',0:.1:.4,...
    'xtick',.8:.2:1.6,'xticklabel',{'' 1, '' 1.4 ''})
% xlabel('MP mean, P')
% ylabel('spike variance, RG')

set(gcf,'color','white')
scrsz = get(0,'ScreenSize');
fxDim=1200;
fyDim=400;
set(gcf,'Position',[1 scrsz(4)-fyDim fxDim fyDim])