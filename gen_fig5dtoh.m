%% Figure 5D-H

clear pars
i=1;
pars{i}.k=20*8;
pars{i}.muP=1.42/4.42;
pars{i}.vP=.2/21;
pars{i}.cP=.95;
pars{i}.m=1.4;

i=2;
pars{i}.k=20*4;
pars{i}.muP=1.42/2.69;
pars{i}.vP=.2/7.9;
pars{i}.cP=.95;
pars{i}.m=1.4;

i=3;
pars{i}.k=20*2;
pars{i}.muP=1.42/1.64;
pars{i}.vP=.2/2.95;
pars{i}.cP=.95;
pars{i}.m=1.4;

i=4;
pars{i}.k=20;
pars{i}.muP=1.42;
pars{i}.vP=.2;
pars{i}.cP=.95;
pars{i}.m=1.4;

i=5;
pars{i}.k=20/2;
pars{i}.muP=1.42*1.64;
pars{i}.vP=.2*2.5;
pars{i}.cP=.95;
pars{i}.m=1.4;

i=6;
pars{i}.k=20/4;
pars{i}.muP=1.42*2.7;
pars{i}.vP=.2*6.4;
pars{i}.cP=.95;
pars{i}.m=1.4;

i=7;
pars{i}.k=20/8;
pars{i}.muP=1.42*4.45;
pars{i}.vP=.2*16;
pars{i}.cP=.95;
pars{i}.m=1.4;

i=8;
pars{i}.k=20/16;
pars{i}.muP=1.42*7.26;
pars{i}.vP=.2*47;
pars{i}.cP=.95;
pars{i}.m=1.4;

i=9;
pars{i}.k=20/32;
pars{i}.muP=1.42*11.9;
pars{i}.vP=.2*132;
pars{i}.cP=.95;
pars{i}.m=1.4;

trialNo=200000;
trialSampleNo=5;

for i=1:9,
    [sP] = gen_spikes_nonlin_poiss([trialNo trialSampleNo], ones(2,1)*pars{i}.muP, [1 pars{i}.cP; pars{i}.cP 1]*pars{i}.vP, pars{i}.k/50, 0, pars{i}.m, 1); 
    cc=corrcoef(sP');
    aa=mean(sP,2);
    mPscK(i)=aa(1);
    aa=var(sP,[],2);
    vPscK(i)=aa(1);
    cPscK(i)=cc(1,2);

    sprintf('mu=%1.2f, var=%1.2f, ff=%1.2f, c=%1.3f\n',mPscK(i),vPscK(i),vPscK(i)/mPscK(i),cPsc(i))
end

%%

figID=21;
figure(figID)
clf

kRg=k./(2.^(-3:5));

mpm=zeros(1,length(pars));
mpv=zeros(1,length(pars));
for i=1:length(pars),
    mpm(i)=pars{i}.muP;
    mpv(i)=pars{i}.vP;
end

subplot(1,5,1)
semilogx(kRg,mpm,'o-','linewidth',2,'color',ones(1,3)*0)
set(gca,'ylim',[0 20],'ytick',0:10:20,'box','off','fontsize',16,...
    'xlim',[min(kRg)/2 max(kRg)*2],'xtick',kRg(end:-2:1),'xticklabel',{kRg(end) '' kRg(5) '' kRg(1)})
ylabel('\mu_{MP}')
xlabel('k')

subplot(1,5,2)
semilogx(kRg,mpv,'o-','linewidth',2,'color',ones(1,3)*0)
set(gca,'ylim',[0 30],'ytick',0:10:30,'box','off','fontsize',16,...
    'xlim',[min(kRg)/2 max(kRg)*2],'xtick',kRg(end:-2:1),'xticklabel',{kRg(end) '' kRg(5) '' kRg(1)})
ylabel('\sigma^2_{MP}')
xlabel('k')

subplot(1,5,3)
semilogx(kRg,mPscK,'o-','linewidth',2,'color',ones(1,3)*0)
set(gca,'ylim',[0 4],'ytick',0:4,'yticklabel',{0 10 20 30 40},'box','off','fontsize',16,...
    'xlim',[min(kRg)/2 max(kRg)*2],'xtick',kRg(end:-2:1),'xticklabel',{kRg(end) '' kRg(5) '' kRg(1)})
ylabel('spike count')
xlabel('k')

subplot(1,5,4)
semilogx(k./(2.^(-3:5)),vPscK./mPscK,'o-','linewidth',2,'color',ones(1,3)*0)
set(gca,'ylim',[0 1.5],'fontsize',16,'box','off',...
    'xlim',[min(kRg)/2 max(kRg)*2],'xtick',kRg(end:-2:1),'xticklabel',{kRg(end) '' kRg(5) '' kRg(1)})
ylabel('Fano factor')
xlabel('k')

subplot(1,5,5)
semilogx(k./(2.^(-3:5)),cPscK,'o-','linewidth',2,'color',ones(1,3)*0)
set(gca,'ylim',[0 .5],'ytick',0:.1:.5,'yticklabel',{0 '' '' '' '' 0.5},'fontsize',16,'box','off',...
    'xlim',[min(kRg)/2 max(kRg)*2],'xtick',kRg(end:-2:1),'xticklabel',{kRg(end) '' kRg(5) '' kRg(1)})
ylabel('spike count correlation')
xlabel('k')

set(gcf,'color','white')
scrsz = get(0,'ScreenSize');
fxDim=1000;
fyDim=170;
set(gcf,'Position',[1 scrsz(4)-fyDim fxDim fyDim])