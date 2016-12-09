%% Figure 5I,J

rg=-1:2/10:2;

iRg=1:9;

mPscRg=zeros(numel(iRg),numel(rg));
vPscRg=zeros(numel(iRg),numel(rg));
cPscRg=zeros(numel(iRg),numel(rg));
for i=1:numel(iRg),
    i
    for n=1:numel(rg),
        [sP] = gen_spikes_nonlin_poiss([trialNo trialSampleNo], ...
            ones(2,1)*pars{i}.muP*rg(n), ...
            [1 pars{i}.cP; pars{i}.cP 1]*pars{i}.vP, ...
            pars{i}.k/50, 0, pars{i}.m, 1); 
        cc=corrcoef(sP');
        aa=mean(sP,2);
        mPscRg(i,n)=aa(1);
        aa=var(sP,[],2);
        vPscRg(i,n)=aa(1);
        cPscRg(i,n)=cc(1,2);
    end
end

%%

figID=23;
figure(figID);
clf

[maxcorr maxcorri]=max(mean(cPscRg,1));
ii=(find(rg==1));
cPscRgm = mean(cPscRg,1);

xl=0.2;

subplot(1,2,1)
hold on
for i=1:numel(iRg),
    plot(rg,mPscRg(i,:),'color',(10-i)*.1*ones(1,3),'linewidth',1)
end
plot(rg(maxcorri),mPscRg(1,maxcorri),'d','color',[.9 0 .1],'markersize',15,'linewidth',3)
plot(rg(ii),mPscRg(1,ii),'o','color',zeros(1,3),'markersize',15,'linewidth',3)
set(gca,'ylim',[0 10],'fontsize',16,'xlim',[min(rg) max(rg)],'ytick',[0:5:10],'yticklabel',{0 50 100},...
    'xtick',-1:2)
ylabel('firing rate')
xlabel('normalized membrane potential')


subplot(1,2,2)
hold on
for i=1:numel(iRg),
    plot(rg,cPscRg(i,:),'color',(10-i)*.1*ones(1,3),'linewidth',1)
end
plot(rg(ii),cPscRgm(1,ii),'o','color',zeros(1,3),'markersize',15,'linewidth',3)
plot(rg(maxcorri),maxcorr,'d','color',[.9 0 .1],'markersize',15,'linewidth',3)
set(gca,'ylim',[0 .5],'fontsize',16,'xlim',[min(rg) max(rg)],'ytick',[0:.25:.5],...
    'xtick',-1:2)
ylabel('spike count correlation')
xlabel('normalized membrane potential')

set(gcf,'color','white')
scrsz = get(0,'ScreenSize');
fxDim=1000;
fyDim=350;
set(gcf,'Position',[1 scrsz(4)-fyDim fxDim fyDim])
