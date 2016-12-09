%% Figure 5A

muRGorig=0.9;
muPorig=1.5;
cRG=.13;
cP=.95;
vRG=2.75;
vP=.2;
k=20;
m=1.4;
trialSampleNo=5;

trialNo=200000;

rg=-2:.2:3;
rg2=-.95:0.1:0.95;

rgl=numel(rg);
rgl2=numel(rg2);

cPscs=zeros(rgl2,rgl);
cRGscs=zeros(rgl2,rgl);
muRGrg=zeros(1,rgl);
muPrg=zeros(1,rgl);
for j=1:rgl2,
    cP=rg2(j);
    for i=1:rgl,
        muRG = muRGorig -.9 +rg(i);
        muP = muPorig -.55 +.6*rg(i);
        muRGrg(i)=muRG;
        muPrg(i)=muP;

        [sP] = gen_spikes_nonlin_poiss([trialNo trialSampleNo], ones(2,1)*muP, [1 cP; cP 1]*vP, k/50, 0, m, 1); 
        cc=corrcoef(sP');
        cPscs(j,i)=cc(1,2);
        mPscs(:,i)=mean(sP,2);
        vPscs(:,i)=var(sP,[],2);
    end
end

%%

figID = 6;
figure(figID)
clf

colormap('gray')

hold on;
for i=1:18
    plot(muPrg,cPscs(i,:),'color',ones(1,3)*(i-1)/30)
end
set(gca,'ylim',[-.2 .2],'fontsize',16,'box','off','xlim',[-.3,2.8],'xtick',0:.5:2.5,'xticklabel',{0 '' 1 '' 2 ''},'ytick',-.2:.1:.2,'yticklabel',{-0.2 '' 0 '' 0.2})
xlabel('membrane potential')
ylabel('spike count correlation')

colorbar('ylim',[0 .6],'ytick',[0 0.3 .6],'yticklabel',{'-0.85' '0' '0.85'})

set(gcf,'color','white')
scrsz = get(0,'ScreenSize');
fxDim=500;
fyDim=400;
set(gcf,'Position',[1 scrsz(4)-fyDim fxDim fyDim])