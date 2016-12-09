%% Figure 5B

muRGorig=0.9;
muPorig=1.5;
cRG=.13;
cP=.95;
vRG=2.75;
vP=.2;
k=20;
m=1.2;
trialSampleNo=5;

trialNo=200000;

rg=-2:.2:3;
rg2=-.4:0.05:0.4;

rgl=numel(rg);
rgl2=numel(rg2);

cPscs=zeros(rgl2,rgl);
cRGscs=zeros(rgl2,rgl);
muRGrg=zeros(1,rgl);
muPrg=zeros(1,rgl);
for j=1:rgl2,
    j
    cRG=rg2(j);
    for i=1:rgl,
        muRG = muRGorig -.9 +rg(i);
        muP = muPorig -.55 +.6*rg(i);
        muRGrg(i)=muRG;
        muPrg(i)=muP;

        [sRG] = gen_spikes_nonlin_nonpoiss([trialNo trialSampleNo], ones(2,1)*muRG, [1 cRG; cRG 1]*vRG, k/50, 0, m, 1);
        cc=corrcoef(sRG');
        cRGscs(j,i)=cc(1,2);
        mRGscs(:,i)=mean(sRG,2);
        vRGscs(:,i)=var(sRG,[],2);
    end
end

%%

figID = 7;
figure(figID)
clf

colormap('gray')

hold on;
for i=1:17
    plot(muRGrg,cRGscs(i,:),'color',ones(1,3)*(i-1)/30)
end
set(gca,'ylim',[-.45 .45],'fontsize',16,'box','off','xlim',[-2.1,3.1],'xtick',-2:1:3,'xticklabel',{-2 '' 0 '' 2 ''},'ytick',-.4:.2:.4,'yticklabel',{-0.4 '' 0 '' 0.4})
xlabel('membrane potential')
ylabel('spike count correlation')

colorbar('ylim',[0 .6],'ytick',[0 0.3 .6],'yticklabel',{'-0.4' '0' '0.4'})

set(gcf,'color','white')
scrsz = get(0,'ScreenSize');
fxDim=500;
fyDim=400;
set(gcf,'Position',[1 scrsz(4)-fyDim fxDim fyDim])