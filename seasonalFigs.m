% Seasonal cycle plots for globe and 3 regions


%% all 12 cases, global seasonal cycle
load('yr10JN100m.mat','JNall*')
load('yr10N100m.mat','Nall*')
load('globalLatlonbasin.mat')
load('stratificationNutrientQ.mat', 'taream')
logN=(lat>0)&(basin~=0);%northern and southern hemisphere
logS=(lat<0)&(basin~=0);
arean=sum(taream(logN));
areas=sum(taream(logS));
areat=sum(taream(basin~=0));
%
prod=squeeze(nansum(-JNall2000*(10*86400*365*1e-3*14*117/16).*repmat(taream,[1 1 10 12 12]),3));%xyz time case
prod1=squeeze(nansum(-JNall2100*(10*86400*365*1e-3*14*117/16).*repmat(taream,[1 1 10 12 12]),3));

for i=1:12
    for j=1:12
        hold1=prod(:,:,i,j);
        prodAn(i,j)=nansum(hold1(logN));
        prodAs(i,j)=nansum(hold1(logS)); 
        hold1=prod1(:,:,i,j);
        prodAn1(i,j)=nansum(hold1(logN));
        prodAs1(i,j)=nansum(hold1(logS)); 
    end
end
prodG0=((prodAn([12 1:11],:))+(prodAs([6:12 1:5],:)))./(areat);
prodG1=((prodAn1([12 1:11],:))+(prodAs1([6:12 1:5],:)))./(areat);
%
figure; 
for i=1:12
    subplot(4,3,i)
    plot(1:12,prodG0(:,i),1:12,prodG1(:,i),'--')
    xlim([0.5 12.5]); ylim([0 44]); %title(num2str(orderparam12(i)))
end
% plot ~dQL
figure;
for i=1:12
    subplot(4,3,i)
    plot(1:12,(-prodG0(:,i)+prodG1(:,i))/mu(i))
    xlim([0.5 12.5]); %title(num2str(orderparam12(i)))
end
%% arctic
arctic=(lat>66.5)&(basin~=0);
areat=sum(taream(arctic));
for i=1:12
    for j=1:12
        hold1=prod(:,:,i,j);
        prodA(i,j)=nansum(hold1(arctic));
        hold1=prod1(:,:,i,j);
        prodA1(i,j)=nansum(hold1(arctic));
    end
end
prodArc0=prodA([12 1:11],:)./(areat);
prodArc1=prodA1([12 1:11],:)./(areat);
figure; 
for i=1:12
    subplot(4,3,i)
    plot(1:12,prodArc0(:,i),1:12,prodArc1(:,i),'--')
    xlim([0.5 12.5]); % title(num2str(orderparam12(i)))
end
% plot ~dQL
figure;
for i=1:12
    subplot(4,3,i)
    plot(1:12,(-prodArc0(:,i)+prodArc1(:,i))/mu(i))
    xlim([0.5 12.5]); %title(num2str(orderparam12(i)))
end
%%
%load NPbudgets
%entrainment1, production1, export1
load('npFluxes100m2100tau30w10.332.mat', 'entrainment1','export1','production1','taream')
entrainb1=entrainment1(:,:,end-35:end);
exportb1=export1(:,:,end-35:end);
prodb1=production1(:,:,end-35:end);%in N units, gN/yr; multiply by 117/106 for C
load('npFluxes100m2000tau30w10.332.mat', 'entrainment1','export1','production1')
entrainb=entrainment1(:,:,end-35:end);
exportb=export1(:,:,end-35:end);
prodb=production1(:,:,end-35:end);
load('npFluxes100m2000tau1yw5.021.mat', 'entrainment1','export1','production1')
entraina=entrainment1(:,:,end-35:end);
exporta=export1(:,:,end-35:end);
proda=production1(:,:,end-35:end);
load('npFluxes100m2100tau1yw5.021.mat', 'entrainment1','export1','production1')
entraina1=entrainment1(:,:,end-35:end);
exporta1=export1(:,:,end-35:end);
proda1=production1(:,:,end-35:end);

%load XMXL %just for PAP/osmosis
XMXL=0.01*(ncread('g.e21.G.T62_g17.param2000.121.pop.h.XMXL.002101-003012.nc','XMXL',[1 1 109],[Inf Inf Inf]));%0.01 to go cm-->m
XMXL1=0.01*(ncread('g.e21.G1850ECO.T62_g17.param2100.121.pop.h.XMXL.002101-003012.nc','XMXL',[1 1 109],[Inf Inf Inf]));
dXMXL=XMXL1-XMXL;

%load stratification
PD=1000*(ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','PD',[1 1 1 109],[Inf Inf 16 Inf])-1);
PD1=1000*(ncread('g.e21.G1850ECO.T62_g17.param2100.123.pop.h.PD.002101-003012.nc','PD',[1 1 1 109],[Inf Inf 16 Inf])-1);
z=ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','z_t')./100;
%strat=squeeze(PD(:,:,16,end-11:end)-PD(:,:,6,end-11:end))./100;%year 10  potential density at 155m-55m/100m
%strat1=squeeze(PD1(:,:,16,end-11:end)-PD1(:,:,6,end-11:end))./100;
%dstrat=strat1-strat;
strat=squeeze(PD(:,:,16,end-11:end)-PD(:,:,1,end-11:end))./-150;  %pd 155m-5m/150
strat1=squeeze(PD1(:,:,16,end-11:end)-PD(:,:,1,end-11:end))./-150; 
clear PD PD1

%load light
load('climateAndRegions.mat', 'light','light1','mld','mld1')
light=light(:,:,end-35:end);
light1=light1(:,:,end-35:end);
mld=mld(:,:,end-35:end);
mld1=mld1(:,:,end-35:end);

%compute Q,L
load('yr10JN100m.mat','JNall*')
load('yr10N100m.mat','Nall*')
load('globalLatlonbasin.mat')
JNall2000=JNall2000(:,:,:,:,[4 11]);
JNall2100=JNall2100(:,:,:,:,[4 11]);
Nall2000=Nall2000(:,:,:,:,[4 11]);
Nall2100=Nall2100(:,:,:,:,[4 11]);

kn(1,1,1,1,1:2)=[0.25 1];
mu(1,1,1,1,1:2)=[0.5 2]./86400;

Q0=Nall2000./(Nall2000+repmat(kn,[320 384 10 12 1]));
Q1=Nall2100./(Nall2100+repmat(kn,[320 384 10 12 1]));
L0=-JNall2000./(repmat(mu,[320 384 10 12 1]).*Q0);
L1=-JNall2100./(repmat(mu,[320 384 10 12 1]).*Q1);

meanQ0=squeeze(nanmean(Q0,3));%depth average
meanQ1=squeeze(nanmean(Q1,3));
meanL0=squeeze(nanmean(L0,3));
meanL1=squeeze(nanmean(L1,3));
%clear Q0 Q1 L0 L1
dQLs=squeeze(nanmean(Q1.*L1-Q0.*L0,3));
LdQs=squeeze(nanmean(L0.*(Q1-Q0),3));
QdLs=squeeze(nanmean(Q0.*(L1-L0),3));
dQdLs=squeeze(nanmean((L1-L0).*(Q1-Q0),3));

QL0=squeeze(nanmean(Q0.*L0,[3 4]));
%% QL decomp regional averages
load('climateAndRegions.mat', 'subtropSPac','osmosis','arctic')
logN=(lat>0)&(basin~=0);%northern and southern hemisphere
logS=(lat<0)&(basin~=0);
fracN=sum(taream(logN))./sum(taream(basin~=0));

% !!need to weight by area!!
areaSSP=sum(taream(subtropSPac));
areaArctic=sum(taream(arctic));
areaOs=sum(taream(osmosis));
areaN=sum(taream(logN));
areaS=sum(taream(logS));

for monthi=1:12
    for k=1:2
    hold1=dQLs(:,:,monthi,k);
    dQLr(monthi,k,1)= nansum(hold1(subtropSPac).*taream(subtropSPac))./areaSSP;
    dQLr(monthi,k,2)= nansum(hold1(arctic).*taream(arctic))./areaArctic;
    dQLr(monthi,k,3)= nansum(hold1(osmosis).*taream(osmosis))./areaOs;
    dQLr(monthi,k,4)= fracN*nansum(hold1(logN).*taream(logN))./areaN;
    dQLr(monthi,k,5)= (1-fracN)*nansum(hold1(logS).*taream(logS))./areaS;
        hold1=LdQs(:,:,monthi,k);
    LdQr(monthi,k,1)= nansum(hold1(subtropSPac).*taream(subtropSPac))./areaSSP;
    LdQr(monthi,k,2)= nansum(hold1(arctic).*taream(arctic))./areaArctic;
    LdQr(monthi,k,3)= nansum(hold1(osmosis).*taream(osmosis))./areaOs;
    LdQr(monthi,k,4)= fracN*nansum(hold1(logN).*taream(logN))./areaN;
    LdQr(monthi,k,5)= (1-fracN)*nansum(hold1(logS).*taream(logS))./areaS;
        hold1=QdLs(:,:,monthi,k);
    QdLr(monthi,k,1)= nansum(hold1(subtropSPac).*taream(subtropSPac))./areaSSP;
    QdLr(monthi,k,2)= nansum(hold1(arctic).*taream(arctic))./areaArctic;
    QdLr(monthi,k,3)= nansum(hold1(osmosis).*taream(osmosis))./areaOs;
    QdLr(monthi,k,4)= fracN*nansum(hold1(logN).*taream(logN))./areaN;
    QdLr(monthi,k,5)= (1-fracN)*nansum(hold1(logS).*taream(logS))./areaS;
        hold1=dQdLs(:,:,monthi,k);
    dQdLr(monthi,k,1)= nansum(hold1(subtropSPac).*taream(subtropSPac))./areaSSP;
    dQdLr(monthi,k,2)= nansum(hold1(arctic).*taream(arctic))./areaArctic;
    dQdLr(monthi,k,3)= nansum(hold1(osmosis).*taream(osmosis))./areaOs;
    dQdLr(monthi,k,4)= fracN*nansum(hold1(logN).*taream(logN))./areaN;
    dQdLr(monthi,k,5)= (1-fracN)*nansum(hold1(logS).*taream(logS))./areaS;
    end
end

%% global 6-month offset cycle computation

logN=(lat>0)&(basin~=0);%northern and southern hemisphere
logS=(lat<0)&(basin~=0);
arean=sum(taream(logN));
areas=sum(taream(logS));
areat=sum(taream(basin~=0));

% !! need to area-weight these too- for each location, not just hemisphere!!

for i=1:36
    hold1=entraina(:,:,i);
    entrainAn(i)=nansum(hold1(logN));
    entrainAs(i)=nansum(hold1(logS));
    hold1=entraina1(:,:,i);
    entrainAn1(i)=nansum(hold1(logN));
    entrainAs1(i)=nansum(hold1(logS));
    hold1=entrainb1(:,:,i);
    entrainBn1(i)=nansum(hold1(logN));
    entrainBs1(i)=nansum(hold1(logS));
    hold1=entrainb(:,:,i);
    entrainBn(i)=nansum(hold1(logN));
    entrainBs(i)=nansum(hold1(logS));

    hold1=exporta(:,:,i);
    exportAn(i)=nansum(hold1(logN));
    exportAs(i)=nansum(hold1(logS));
    hold1=exporta1(:,:,i);
    exportAn1(i)=nansum(hold1(logN));
    exportAs1(i)=nansum(hold1(logS));
    hold1=exportb1(:,:,i);
    exportBn1(i)=nansum(hold1(logN));
    exportBs1(i)=nansum(hold1(logS));
    hold1=exportb(:,:,i);
    exportBn(i)=nansum(hold1(logN));
    exportBs(i)=nansum(hold1(logS));
    
    hold1=proda(:,:,i);
    prodAn(i)=nansum(hold1(logN));
    prodAs(i)=nansum(hold1(logS));
    hold1=proda1(:,:,i);
    prodAn1(i)=nansum(hold1(logN));
    prodAs1(i)=nansum(hold1(logS));
    hold1=prodb1(:,:,i);
    prodBn1(i)=nansum(hold1(logN));
    prodBs1(i)=nansum(hold1(logS));
    hold1=prodb(:,:,i);
    prodBn(i)=nansum(hold1(logN));
    prodBs(i)=nansum(hold1(logS));
    
    hold1=mld(:,:,i);
    mldn(i)=nansum(hold1(logN).*taream(logN));
    mlds(i)=nansum(hold1(logS).*taream(logS));
    hold1=mld1(:,:,i);
    mldn1(i)=nansum(hold1(logN).*taream(logN));
    mlds1(i)=nansum(hold1(logS).*taream(logS));
    
    hold1=light(:,:,i);
    lightn(i)=nansum(hold1(logN).*taream(logN));
    lights(i)=nansum(hold1(logS).*taream(logS));
    hold1=light1(:,:,i);
    lightn1(i)=nansum(hold1(logN).*taream(logN));
    lights1(i)=nansum(hold1(logS).*taream(logS));
     
    if i<13
        hold1=strat(:,:,i);
        stratn(i)=nansum(hold1(logN).*taream(logN));
        strats(i)=nansum(hold1(logS).*taream(logS));
        hold1=strat1(:,:,i);
        stratn1(i)=nansum(hold1(logN).*taream(logN));
        strats1(i)=nansum(hold1(logS).*taream(logS));
    
        hold1=meanQ0(:,:,i,1);
        QAn(i)=nansum(hold1(logN).*taream(logN));
        QAs(i)=nansum(hold1(logS).*taream(logS));
        hold1=meanQ0(:,:,i,2);
        QBn(i)=nansum(hold1(logN).*taream(logN));
        QBs(i)=nansum(hold1(logS).*taream(logS));
        hold1=meanQ1(:,:,i,1);
        QAn1(i)=nansum(hold1(logN).*taream(logN));
        QAs1(i)=nansum(hold1(logS).*taream(logS));
        hold1=meanQ1(:,:,i,2);
        QBn1(i)=nansum(hold1(logN).*taream(logN));
        QBs1(i)=nansum(hold1(logS).*taream(logS));
        
        hold1=meanL0(:,:,i,1);
        LAn(i)=nansum(hold1(logN).*taream(logN));
        LAs(i)=nansum(hold1(logS).*taream(logS));
        hold1=meanL0(:,:,i,2);
        LBn(i)=nansum(hold1(logN).*taream(logN));
        LBs(i)=nansum(hold1(logS).*taream(logS));
        hold1=meanL1(:,:,i,1);
        LAn1(i)=nansum(hold1(logN).*taream(logN));
        LAs1(i)=nansum(hold1(logS).*taream(logS));
        hold1=meanL1(:,:,i,2);
        LBn1(i)=nansum(hold1(logN).*taream(logN));
        LBs1(i)=nansum(hold1(logS).*taream(logS));
    end
end

%global seasonal cycle
entrainGb0=((entrainBn([36 25:35])+entrainBn([24 13:23])+entrainBn([12 1:11]))+(entrainBs([30:36 25:29])+entrainBs([18:24 13:17])+entrainBs([6:12 1:5])))*0.365*14*117/(16*3*areat);
entrainGb1=((entrainBn1([36 25:35])+entrainBn1([24 13:23])+entrainBn1([12 1:11]))+(entrainBs1([30:36 25:29])+entrainBs1([18:24 13:17])+entrainBs1([6:12 1:5])))*0.365*14*117/(16*3*areat);
entrainGa0=((entrainAn([36 25:35])+entrainAn([24 13:23])+entrainAn([12 1:11]))+(entrainAs([30:36 25:29])+entrainAs([18:24 13:17])+entrainAs([6:12 1:5])))*0.365*14*117/(16*3*areat);
entrainGa1=((entrainAn1([36 25:35])+entrainAn1([24 13:23])+entrainAn1([12 1:11]))+(entrainAs1([30:36 25:29])+entrainAs1([18:24 13:17])+entrainAs1([6:12 1:5])))*0.365*14*117/(16*3*areat);

exportGb0=((exportBn([36 25:35])+exportBn([24 13:23])+exportBn([12 1:11]))+(exportBs([30:36 25:29])+exportBs([18:24 13:17])+exportBs([6:12 1:5])))*0.365*14*117/(16*3*areat);
exportGb1=((exportBn1([36 25:35])+exportBn1([24 13:23])+exportBn1([12 1:11]))+(exportBs1([30:36 25:29])+exportBs1([18:24 13:17])+exportBs1([6:12 1:5])))*0.365*14*117/(16*3*areat);
exportGa0=((exportAn([36 25:35])+exportAn([24 13:23])+exportAn([12 1:11]))+(exportAs([30:36 25:29])+exportAs([18:24 13:17])+exportAs([6:12 1:5])))*0.365*14*117/(16*3*areat);
exportGa1=((exportAn1([36 25:35])+exportAn1([24 13:23])+exportAn1([12 1:11]))+(exportAs1([30:36 25:29])+exportAs1([18:24 13:17])+exportAs1([6:12 1:5])))*0.365*14*117/(16*3*areat);

prodGb0=((prodBn([36 25:35])+prodBn([24 13:23])+prodBn([12 1:11]))+(prodBs([30:36 25:29])+prodBs([18:24 13:17])+prodBs([6:12 1:5])))*0.365*14*117/(16*3*areat);
prodGb1=((prodBn1([36 25:35])+prodBn1([24 13:23])+prodBn1([12 1:11]))+(prodBs1([30:36 25:29])+prodBs1([18:24 13:17])+prodBs1([6:12 1:5])))*0.365*14*117/(16*3*areat);
prodGa0=((prodAn([36 25:35])+prodAn([24 13:23])+prodAn([12 1:11]))+(prodAs([30:36 25:29])+prodAs([18:24 13:17])+prodAs([6:12 1:5])))*0.365*14*117/(16*3*areat);
prodGa1=((prodAn1([36 25:35])+prodAn1([24 13:23])+prodAn1([12 1:11]))+(prodAs1([30:36 25:29])+prodAs1([18:24 13:17])+prodAs1([6:12 1:5])))*0.365*14*117/(16*3*areat);


stratG1=((stratn1([12 1:11]))+(strats1([12 1:11])))/(3*areat);
stratG0=((stratn([12 1:11]))+(strats([12 1:11])))/(3*areat);
mldG1=((mldn1([36 25:35])+mldn1([24 13:23])+mldn1([12 1:11]))+(mlds1([30:36 25:29])+mlds1([18:24 13:17])+mlds1([6:12 1:5])))/(3*areat);
mldG0=((mldn([36 25:35])+mldn([24 13:23])+mldn([12 1:11]))+(mlds([30:36 25:29])+mlds([18:24 13:17])+mlds([6:12 1:5])))/(3*areat);

lightG1=((lightn1([36 25:35])+lightn1([24 13:23])+lightn1([12 1:11]))+(lights1([30:36 25:29])+lights1([18:24 13:17])+lights1([6:12 1:5])))/(3*areat);
lightG0=((lightn([36 25:35])+lightn([24 13:23])+lightn([12 1:11]))+(lights([30:36 25:29])+lights([18:24 13:17])+lights([6:12 1:5])))/(3*areat);


QGa0=(QAn([12 1:11])+QAs([6:12 1:5]))/areat;
QGb0=(QBn([12 1:11])+QBs([6:12 1:5]))/areat;
QGa1=(QAn1([12 1:11])+QAs1([6:12 1:5]))/areat;
QGb1=(QBn1([12 1:11])+QBs1([6:12 1:5]))/areat;

LGa0=(LAn([12 1:11])+LAs([6:12 1:5]))/areat;
LGb0=(LBn([12 1:11])+LBs([6:12 1:5]))/areat;
LGa1=(LAn1([12 1:11])+LAs1([6:12 1:5]))/areat;
LGb1=(LBn1([12 1:11])+LBs1([6:12 1:5]))/areat;

dqlGa=dQLr([12 1:11],1,4)+dQLr([6:12 1:5],1,5);
qdlGa=QdLr([12 1:11],1,4)+QdLr([6:12 1:5],1,5);
ldqGa=LdQr([12 1:11],1,4)+LdQr([6:12 1:5],1,5);
dqdlGa=dQdLr([12 1:11],1,4)+dQdLr([6:12 1:5],1,5);
dqlGb=dQLr([12 1:11],2,4)+dQLr([6:12 1:5],2,5);
qdlGb=QdLr([12 1:11],2,4)+QdLr([6:12 1:5],2,5);
ldqGb=LdQr([12 1:11],2,4)+LdQr([6:12 1:5],2,5);
dqdlGb=dQdLr([12 1:11],2,4)+dQdLr([6:12 1:5],2,5);

globalmaxQL0=squeeze(max(max(QL0)));
globalmeanQL0=squeeze(areaweightedmean(areaweightedmean(QL0,taream,2),mean(taream,2),1));
%% global plot
%entrain,prod,export; strat,mld; QL decomp; light

figure;
subplot(3,2,1) %entrain
plot(1:12,entrainGa0,1:12,entrainGb0); hold all
plot(1:12,entrainGa1,'--','Color',[0 0.45 0.74])
plot(1:12,entrainGb1,'--','Color',[0.85 0.33 0.1])
set(gca,'XTick',1:12); xlim([1 12]); ylim([0 40])
set(gca,'XTickLabels',{})
ylabel('gC/m^2 y'); xlabel('(a)')
set(gca,'fontsize',12)

subplot(3,2,2) %prod
plot(1:12,prodGa0,1:12,prodGb0); hold all
plot(1:12,prodGa1,'--','Color',[0 0.45 0.74])
plot(1:12,prodGb1,'--','Color',[0.85 0.33 0.1])
set(gca,'XTick',1:12); xlim([1 12]); ylim([0 40])
set(gca,'XTickLabels',{});%set(gca,'XTickLabel',{'J\newlineJ','F\newlineA','M\newlineS','A\newlineO','M\newlineN','J\newlineD','J\newlineJ','A\newlineF','S\newlineM','O\newlineA','N\newlineM','D\newlineJ'}) 
ylabel('gC/m^2 y'); xlabel('(b)')
set(gca,'fontsize',12)
legend('slow 2000','fast 2000','slow 2100','fast 2100')

% subplot(3,2,5) %export
% plot(1:12,exportGa0,1:12,exportGb0); hold all
% plot(1:12,exportGa1,'--','Color',[0 0.45 0.74])
% plot(1:12,exportGb1,'--','Color',[0.85 0.33 0.1])
% set(gca,'XTick',1:12); xlim([1 12]); ylim([-40 0])
% set(gca,'XTickLabel',{'J\newlineJ','F\newlineA','M\newlineS','A\newlineO','M\newlineN','J\newlineD','J\newlineJ','A\newlineF','S\newlineM','O\newlineA','N\newlineM','D\newlineJ'}) 
% ylabel('gC/m^2 y'); xlabel('(c)')
% set(gca,'fontsize',12); 

subplot(3,2,5) %mld, no strat
%yyaxis 'left'
% set(gca,'YColor','k');
plot(1:12,mldG0/100,'k'); hold on; plot(1:12,mldG1/100,'k--')
ylabel('MLD')
%yyaxis('right')
%set(gca,'YColor','r');
%plot(1:12,stratG0*150,'r'); hold on
%plot(1:12,stratG1*150,'r--')
%ylabel('\Delta \sigma_\theta')
set(gca,'XTick',1:12); xlim([1 12])
set(gca,'XTickLabel',{'J\newlineJ','F\newlineA','M\newlineS','A\newlineO','M\newlineN','J\newlineD','J\newlineJ','A\newlineF','S\newlineM','O\newlineA','N\newlineM','D\newlineJ'}) 
 xlabel('(e)')
 set(gca,'fontsize',12); legend('2000','2100')

subplot(3,2,3) %QL decomp
plot(1:12,dqlGa/globalmaxQL0(1),'k','Linewidth',3); hold on
plot(1:12,ldqGa/globalmaxQL0(1),'r'); plot(1:12,qdlGa/globalmaxQL0(1),'c'); plot(1:12,dqdlGa/globalmaxQL0(1),'Color',[0.49 0.18 0.56]);
set(gca,'XTick',1:12); xlim([1 12]);set(gca,'XTickLabels',{})
xlabel('(c)')
set(gca,'fontsize',12)

subplot(3,2,4) %Ql decomp
plot(1:12,dqlGb/globalmaxQL0(2),'k','Linewidth',3); hold on
plot(1:12,ldqGb/globalmaxQL0(2),'r'); plot(1:12,qdlGb/globalmaxQL0(2),'c'); plot(1:12,dqdlGb/globalmaxQL0(2),'Color',[0.49 0.18 0.56]);
set(gca,'XTick',1:12); xlim([1 12])
%set(gca,'XTickLabels',{})
set(gca,'XTickLabel',{'J\newlineJ','F\newlineA','M\newlineS','A\newlineO','M\newlineN','J\newlineD','J\newlineJ','A\newlineF','S\newlineM','O\newlineA','N\newlineM','D\newlineJ'}) 
xlabel('(d)'); legend('\Delta QL','L \DeltaQ','Q \DeltaL','\DeltaQ \DeltaL')
set(gca,'fontsize',12)


%% subtropical south pacific computation
logSSP=subtropSPac;
areat=sum(taream(logSSP));


for i=1:36
    hold1=entraina(:,:,i);
    entrainAn(i)=nansum(hold1(logSSP))./areat;
    hold1=entraina1(:,:,i);
    entrainAn1(i)=nansum(hold1(logSSP))./areat;
    hold1=entrainb1(:,:,i);
    entrainBn1(i)=nansum(hold1(logSSP))./areat;
    hold1=entrainb(:,:,i);
    entrainBn(i)=nansum(hold1(logSSP))./areat;

    hold1=exporta(:,:,i);
    exportAn(i)=nansum(hold1(logSSP))./areat;
    hold1=exporta1(:,:,i);
    exportAn1(i)=nansum(hold1(logSSP))./areat;
    hold1=exportb1(:,:,i);
    exportBn1(i)=nansum(hold1(logSSP))./areat;
    hold1=exportb(:,:,i);
    exportBn(i)=nansum(hold1(logSSP))./areat;
    
    hold1=proda(:,:,i);
    prodAn(i)=nansum(hold1(logSSP))./areat;
    hold1=proda1(:,:,i);
    prodAn1(i)=nansum(hold1(logSSP))./areat;
    hold1=prodb1(:,:,i);
    prodBn1(i)=nansum(hold1(logSSP))./areat;
    hold1=prodb(:,:,i);
    prodBn(i)=nansum(hold1(logSSP))./areat;

    if i<13
        hold1=meanQ0(:,:,i,1);
        QAn(i)=nansum(hold1(logSSP).*taream(logSSP));
        hold1=meanQ0(:,:,i,2);
        QBn(i)=nansum(hold1(logSSP).*taream(logSSP));
        hold1=meanQ1(:,:,i,1);
        QAn1(i)=nansum(hold1(logSSP).*taream(logSSP));
        hold1=meanQ1(:,:,i,2);
        QBn1(i)=nansum(hold1(logSSP).*taream(logSSP));
        
        hold1=meanL0(:,:,i,1);
        LAn(i)=nansum(hold1(logSSP).*taream(logSSP));
        hold1=meanL0(:,:,i,2);
        LBn(i)=nansum(hold1(logSSP).*taream(logSSP));
        hold1=meanL1(:,:,i,1);
        LAn1(i)=nansum(hold1(logSSP).*taream(logSSP));
        hold1=meanL1(:,:,i,2);
        LBn1(i)=nansum(hold1(logSSP).*taream(logSSP));
    end
end

entrainSSPb0=(entrainBn([36 25:35])+entrainBn([24 13:23])+entrainBn([12 1:11]))*0.365*14*117/(16*3);
entrainSSPb1=(entrainBn1([36 25:35])+entrainBn1([24 13:23])+entrainBn1([12 1:11]))*0.365*14*117/(16*3);
entrainSSPa0=(entrainAn([36 25:35])+entrainAn([24 13:23])+entrainAn([12 1:11]))*0.365*14*117/(16*3);
entrainSSPa1=(entrainAn1([36 25:35])+entrainAn1([24 13:23])+entrainAn1([12 1:11]))*0.365*14*117/(16*3);

exportSSPb0=(exportBn([36 25:35])+exportBn([24 13:23])+exportBn([12 1:11]))*0.365*14*117/(16*3);
exportSSPb1=(exportBn1([36 25:35])+exportBn1([24 13:23])+exportBn1([12 1:11]))*0.365*14*117/(16*3);
exportSSPa0=(exportAn([36 25:35])+exportAn([24 13:23])+exportAn([12 1:11]))*0.365*14*117/(16*3);
exportSSPa1=(exportAn1([36 25:35])+exportAn1([24 13:23])+exportAn1([12 1:11]))*0.365*14*117/(16*3);

prodSSPb0=(prodBn([36 25:35])+prodBn([24 13:23])+prodBn([12 1:11]))*0.365*14*117/(16*3);
prodSSPb1=(prodBn1([36 25:35])+prodBn1([24 13:23])+prodBn1([12 1:11]))*0.365*14*117/(16*3);
prodSSPa0=(prodAn([36 25:35])+prodAn([24 13:23])+prodAn([12 1:11]))*0.365*14*117/(16*3);
prodSSPa1=(prodAn1([36 25:35])+prodAn1([24 13:23])+prodAn1([12 1:11]))*0.365*14*117/(16*3);

dqlSSPa=dQLr([12 1:11],1,1);
qdlSSPa=QdLr([12 1:11],1,1);
ldqSSPa=LdQr([12 1:11],1,1);
dqdlSSPa=dQdLr([12 1:11],1,1);
dqlSSPb=dQLr([12 1:11],2,1);
qdlSSPb=QdLr([12 1:11],2,1);
ldqSSPb=LdQr([12 1:11],2,1);
dqdlSSPb=dQdLr([12 1:11],2,1);

for i=1:2
    holdvar=QL0(:,:,i);
    QLmaxSSP(i)=max(holdvar(logSSP));
end
%% SSP plot
% Nflux, production, P flux down, growth limit changes 
%ADD entrain change cause??
%switch month order to 7:12 1:6??
figure;
subplot(3,1,1) %entrain 
plot(1:12,entrainSSPa0,1:12,entrainSSPb0); hold all
plot(1:12,entrainSSPa1,'--','Color',[0 0.45 0.74])
plot(1:12,entrainSSPb1,'--','Color',[0.85 0.33 0.1])
set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
set(gca,'XTickLabels',{})
ylabel('gC/m^2 y'); xlabel('(a)')
subplot(3,1,2) %prod
plot(1:12,prodSSPa0,1:12,prodSSPb0); hold all
plot(1:12,prodSSPa1,'--','Color',[0 0.45 0.74])
plot(1:12,prodSSPb1,'--','Color',[0.85 0.33 0.1])
set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
set(gca,'XTickLabels',{})
ylabel('gC/m^2 y'); xlabel('(b)')
% subplot(4,1,3) %export
% plot(1:12,exportSSPa0,1:12,exportSSPb0); hold all
% plot(1:12,exportSSPa1,'--','Color',[0 0.45 0.74])
% plot(1:12,exportSSPb1,'--','Color',[0.85 0.33 0.1])
% set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
% ylabel('gC/m^2 y'); xlabel('(c)')
% set(gca,'XTickLabels',{})
legend('slow 2000','fast 2000','slow 2100','fast 2100')
subplot(3,1,3) %QL decomp
plot(1:12,dqlSSPa/QLmaxSSP(1),1:12,dqlSSPb/QLmaxSSP(2)); hold all;
plot(1:12,ldqSSPa/QLmaxSSP(1),'-o','Color',[0 0.45 0.74])
plot(1:12,ldqSSPb/QLmaxSSP(2),'-o','Color',[0.85 0.33 0.1])
legend('slow \Delta QL','fast \Delta QL','slow L \DeltaQ','fast L \DeltaQ')
set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
ylabel('nondim'); xlabel('(c)')
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'}) 
%%
figure; plot(100*(ldqSSPa-dqlSSPa)./(dqlSSPa)); hold on; plot(100*(ldqSSPb-dqlSSPb)./(dqlSSPb))
%% Arctic computation
% Nflux, production, P flux down, growth limits, incoming light
arctic=(lat>66.5)&(basin~=0);
areat=sum(taream(arctic));
logSS=arctic;

for i=1:36
    hold1=entraina(:,:,i);
    entrainAn(i)=nansum(hold1(logSS))./areat;
    hold1=entraina1(:,:,i);
    entrainAn1(i)=nansum(hold1(logSS))./areat;
    hold1=entrainb1(:,:,i);
    entrainBn1(i)=nansum(hold1(logSS))./areat;
    hold1=entrainb(:,:,i);
    entrainBn(i)=nansum(hold1(logSS))./areat;

    hold1=exporta(:,:,i);
    exportAn(i)=nansum(hold1(logSS))./areat;
    hold1=exporta1(:,:,i);
    exportAn1(i)=nansum(hold1(logSS))./areat;
    hold1=exportb1(:,:,i);
    exportBn1(i)=nansum(hold1(logSS))./areat;
    hold1=exportb(:,:,i);
    exportBn(i)=nansum(hold1(logSS))./areat;
    
    hold1=proda(:,:,i);
    prodAn(i)=nansum(hold1(logSS))./areat;
    hold1=proda1(:,:,i);
    prodAn1(i)=nansum(hold1(logSS))./areat;
    hold1=prodb1(:,:,i);
    prodBn1(i)=nansum(hold1(logSS))./areat;
    hold1=prodb(:,:,i);
    prodBn(i)=nansum(hold1(logSS))./areat;
 
    if i<13
        hold1=meanQ0(:,:,i,1);
        QAn(i)=nansum(hold1(logSS).*taream(logSS))./areat;
        hold1=meanQ0(:,:,i,2);
        QBn(i)=nansum(hold1(logSS).*taream(logSS))./areat;
        hold1=meanQ1(:,:,i,1);
        QAn1(i)=nansum(hold1(logSS).*taream(logSS))./areat;
        hold1=meanQ1(:,:,i,2);
        QBn1(i)=nansum(hold1(logSS).*taream(logSS))./areat;
        
        hold1=meanL0(:,:,i,1);
        LAn(i)=nansum(hold1(logSS).*taream(logSS))./areat;
        hold1=meanL0(:,:,i,2);
        LBn(i)=nansum(hold1(logSS).*taream(logSS))./areat;
        hold1=meanL1(:,:,i,1);
        LAn1(i)=nansum(hold1(logSS).*taream(logSS))./areat;
        hold1=meanL1(:,:,i,2);
        LBn1(i)=nansum(hold1(logSS).*taream(logSS))./areat;
        
        hold1=light(:,:,i);
        lightA(i,2)=nansum(hold1(logSS).*taream(logSS))./areat;
        lightA(i,1)=max(hold1(logSS));
        lightA(i,3)=min(hold1(logSS));
        hold1=light1(:,:,i);
        lightA1(i,2)=nansum(hold1(logSS).*taream(logSS))./areat;
        lightA1(i,1)=max(hold1(logSS));
        lightA1(i,3)=min(hold1(logSS));
        

    end
end

entrainAb0=(entrainBn([36 25:35])+entrainBn([24 13:23])+entrainBn([12 1:11]))*0.365*14*117/(16*3);
entrainAb1=(entrainBn1([36 25:35])+entrainBn1([24 13:23])+entrainBn1([12 1:11]))*0.365*14*117/(16*3);
entrainAa0=(entrainAn([36 25:35])+entrainAn([24 13:23])+entrainAn([12 1:11]))*0.365*14*117/(16*3);
entrainAa1=(entrainAn1([36 25:35])+entrainAn1([24 13:23])+entrainAn1([12 1:11]))*0.365*14*117/(16*3);

exportAb0=(exportBn([36 25:35])+exportBn([24 13:23])+exportBn([12 1:11]))*0.365*14*117/(16*3);
exportAb1=(exportBn1([36 25:35])+exportBn1([24 13:23])+exportBn1([12 1:11]))*0.365*14*117/(16*3);
exportAa0=(exportAn([36 25:35])+exportAn([24 13:23])+exportAn([12 1:11]))*0.365*14*117/(16*3);
exportAa1=(exportAn1([36 25:35])+exportAn1([24 13:23])+exportAn1([12 1:11]))*0.365*14*117/(16*3);

prodAb0=(prodBn([36 25:35])+prodBn([24 13:23])+prodBn([12 1:11]))*0.365*14*117/(16*3);
prodAb1=(prodBn1([36 25:35])+prodBn1([24 13:23])+prodBn1([12 1:11]))*0.365*14*117/(16*3);
prodAa0=(prodAn([36 25:35])+prodAn([24 13:23])+prodAn([12 1:11]))*0.365*14*117/(16*3);
prodAa1=(prodAn1([36 25:35])+prodAn1([24 13:23])+prodAn1([12 1:11]))*0.365*14*117/(16*3);

dqlAa=dQLr([12 1:11],1,2);
qdlAa=QdLr([12 1:11],1,2);
ldqAa=LdQr([12 1:11],1,2);
dqdlAa=dQdLr([12 1:11],1,2);
dqlAb=dQLr([12 1:11],2,2);
qdlAb=QdLr([12 1:11],2,2);
ldqAb=LdQr([12 1:11],2,2);
dqdlAb=dQdLr([12 1:11],2,2);

for i=1:2
    holdvar=QL0(:,:,i);
    QLmaxArc(i)=max(holdvar(logSS));
end
%% arctic plot

figure;
subplot(2,2,1) %entrain
plot(1:12,entrainAa0,1:12,entrainAb0); hold all
plot(1:12,entrainAa1,'--','Color',[0 0.45 0.74])
plot(1:12,entrainAb1,'--','Color',[0.85 0.33 0.1])
set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
set(gca,'XTickLabels',{})
ylabel('gC/m^2 y'); xlabel('(a)')
subplot(2,2,2) %QL slow
plot(1:12,dqlAa/QLmaxArc(1),'k','Linewidth',3); hold on
plot(1:12,ldqAa/QLmaxArc(1),'r'); plot(1:12,qdlAa/QLmaxArc(1),'c'); plot(1:12,dqdlAa/QLmaxArc(1),'Color',[0.49 0.18 0.56]);
set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
set(gca,'XTickLabels',{}); xlabel('(c)')
subplot(2,2,3) %prod
plot(1:12,prodAa0,1:12,prodAb0); hold all
plot(1:12,prodAa1,'--','Color',[0 0.45 0.74])
plot(1:12,prodAb1,'--','Color',[0.85 0.33 0.1])
set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
%set(gca,'XTickLabels',{})
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'}) 
ylabel('gC/m^2 y'); xlabel('(b)')
legend('slow 2000','fast 2000','slow 2100','fast 2100')
subplot(2,2,4) %QL fast
plot(1:12,dqlAb/QLmaxArc(2),'k','Linewidth',3); hold on
plot(1:12,ldqAb/QLmaxArc(2),'r'); plot(1:12,qdlAb/QLmaxArc(2),'c'); plot(1:12,dqdlAb/QLmaxArc(2),'Color',[0.49 0.18 0.56]);
set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
%set(gca,'XTickLabels',{}); 
xlabel('(d)'); legend('\Delta QL','L \DeltaQ','Q \DeltaL','\DeltaQ \DeltaL')
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'}) 
% subplot(3,2,5) %export
% plot(1:12,exportAa0,1:12,exportAb0); hold all
% plot(1:12,exportAa1,'--','Color',[0 0.45 0.74])
% plot(1:12,exportAb1,'--','Color',[0.85 0.33 0.1])
% set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
% set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'}) 
% ylabel('gC/m^2 y'); xlabel('(c)'); legend('slow 2000','fast 2000','slow 2100','fast 2100')
% subplot(3,2,6) %light
% hold on;
% f1=fill([1:12 12:-1:1],[lightA([12 1:11],1).' lightA([11:-1:1 12],3).'],'g','LineStyle','none');
% alpha(f1,0.5)
% f2=fill([1:12 12:-1:1],[lightA1([12 1:11],1).' lightA1([11:-1:1 12],3).'],[0.5 0.5 0.5],'LineStyle','none');
% alpha(f2,0.5)
% plot(1:12,lightA([12 1:11],2),'Color',[0.1 0.6 0.1],'LineWidth',2)
% plot(1:12,lightA1([12 1:11],2),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
%set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
%set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'}) 
%xlabel('(f)'); ylabel('W/m^2'); legend('max 2000','2100')
%% Porcupine Abyssal Plane computation
% Nflux, production, P flux down, growth limit changes, MLD as XMXL
%maskosmosis=(lon>-26)&(lon<-16)&(lat>40)&(lat<50);
logSS=osmosis;
areat=sum(taream(logSS));

for i=1:36
    hold1=entraina(:,:,i);
    entrainAn(i)=nansum(hold1(logSS))./areat;
    hold1=entraina1(:,:,i);
    entrainAn1(i)=nansum(hold1(logSS))./areat;
    hold1=entrainb1(:,:,i);
    entrainBn1(i)=nansum(hold1(logSS))./areat;
    hold1=entrainb(:,:,i);
    entrainBn(i)=nansum(hold1(logSS))./areat;

    hold1=exporta(:,:,i);
    exportAn(i)=nansum(hold1(logSS))./areat;
    hold1=exporta1(:,:,i);
    exportAn1(i)=nansum(hold1(logSS))./areat;
    hold1=exportb1(:,:,i);
    exportBn1(i)=nansum(hold1(logSS))./areat;
    hold1=exportb(:,:,i);
    exportBn(i)=nansum(hold1(logSS))./areat;
    
    hold1=proda(:,:,i);
    prodAn(i)=nansum(hold1(logSS))./areat;
    hold1=proda1(:,:,i);
    prodAn1(i)=nansum(hold1(logSS))./areat;
    hold1=prodb1(:,:,i);
    prodBn1(i)=nansum(hold1(logSS))./areat;
    hold1=prodb(:,:,i);
    prodBn(i)=nansum(hold1(logSS))./areat;
 
    if i<13
        hold1=meanQ0(:,:,i,1);
        QAn(i)=nansum(hold1(logSS).*taream(logSS));
        hold1=meanQ0(:,:,i,2);
        QBn(i)=nansum(hold1(logSS).*taream(logSS));
        hold1=meanQ1(:,:,i,1);
        QAn1(i)=nansum(hold1(logSS).*taream(logSS));
        hold1=meanQ1(:,:,i,2);
        QBn1(i)=nansum(hold1(logSS).*taream(logSS));
        
        hold1=meanL0(:,:,i,1);
        LAn(i)=nansum(hold1(logSS).*taream(logSS));
        hold1=meanL0(:,:,i,2);
        LBn(i)=nansum(hold1(logSS).*taream(logSS));
        hold1=meanL1(:,:,i,1);
        LAn1(i)=nansum(hold1(logSS).*taream(logSS));
        hold1=meanL1(:,:,i,2);
        LBn1(i)=nansum(hold1(logSS).*taream(logSS));
        
        hold1=XMXL(:,:,i);
        xmxlP(i,2)=nansum(hold1(logSS).*taream(logSS))./areat;
        xmxlP(i,1)=max(hold1(logSS));
        xmxlP(i,3)=min(hold1(logSS));
        hold1=XMXL1(:,:,i);
        xmxlP1(i,2)=nansum(hold1(logSS).*taream(logSS))./areat;
        xmxlP1(i,1)=max(hold1(logSS));
        xmxlP1(i,3)=min(hold1(logSS));
    end
end

entrainPb0=(entrainBn([36 25:35])+entrainBn([24 13:23])+entrainBn([12 1:11]))*0.365*14*117/(16*3);
entrainPb1=(entrainBn1([36 25:35])+entrainBn1([24 13:23])+entrainBn1([12 1:11]))*0.365*14*117/(16*3);
entrainPa0=(entrainAn([36 25:35])+entrainAn([24 13:23])+entrainAn([12 1:11]))*0.365*14*117/(16*3);
entrainPa1=(entrainAn1([36 25:35])+entrainAn1([24 13:23])+entrainAn1([12 1:11]))*0.365*14*117/(16*3);

exportPb0=(exportBn([36 25:35])+exportBn([24 13:23])+exportBn([12 1:11]))*0.365*14*117/(16*3);
exportPb1=(exportBn1([36 25:35])+exportBn1([24 13:23])+exportBn1([12 1:11]))*0.365*14*117/(16*3);
exportPa0=(exportAn([36 25:35])+exportAn([24 13:23])+exportAn([12 1:11]))*0.365*14*117/(16*3);
exportPa1=(exportAn1([36 25:35])+exportAn1([24 13:23])+exportAn1([12 1:11]))*0.365*14*117/(16*3);

prodPb0=(prodBn([36 25:35])+prodBn([24 13:23])+prodBn([12 1:11]))*0.365*14*117/(16*3);
prodPb1=(prodBn1([36 25:35])+prodBn1([24 13:23])+prodBn1([12 1:11]))*0.365*14*117/(16*3);
prodPa0=(prodAn([36 25:35])+prodAn([24 13:23])+prodAn([12 1:11]))*0.365*14*117/(16*3);
prodPa1=(prodAn1([36 25:35])+prodAn1([24 13:23])+prodAn1([12 1:11]))*0.365*14*117/(16*3);

dqlPa=dQLr([12 1:11],1,3);
qdlPa=QdLr([12 1:11],1,3);
ldqPa=LdQr([12 1:11],1,3);
dqdlPa=dQdLr([12 1:11],1,3);
dqlPb=dQLr([12 1:11],2,3);
qdlPb=QdLr([12 1:11],2,3);
ldqPb=LdQr([12 1:11],2,3);
dqdlPb=dQdLr([12 1:11],2,3);

for i=1:2
    holdvar=QL0(:,:,i);
    QLmaxPAP(i)=max(holdvar(logSS));
end
%% PAP plot

figure;
subplot(4,1,1) %entrain
plot(1:12,entrainPa0,1:12,entrainPb0); hold all
plot(1:12,entrainPa1,'--','Color',[0 0.45 0.74])
plot(1:12,entrainPb1,'--','Color',[0.85 0.33 0.1])
set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
set(gca,'XTickLabels',{})
ylabel('gC/m^2 y'); xlabel('(a)')

% subplot(3,2,2) %QL slow
% plot(1:12,dqlPa,'Linewidth',3); hold on
% plot(1:12,ldqPa); plot(1:12,qdlPa); plot(1:12,dqdlPa);
% set(gca,'XTick',1:12); xlim([1 12])
% set(gca,'XTickLabels',{}); xlabel('(d)')

subplot(4,1,2) %prod
plot(1:12,prodPa0,1:12,prodPb0); hold all
plot(1:12,prodPa1,'--','Color',[0 0.45 0.74])
plot(1:12,prodPb1,'--','Color',[0.85 0.33 0.1])
set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
set(gca,'XTickLabels',{})
ylabel('gC/m^2 y'); xlabel('(b)'); legend('slow 2000','fast 2000','slow 2100','fast 2100')

% subplot(3,2,4) %QL fast
% plot(1:12,dqlPb,'Linewidth',3); hold on
% plot(1:12,ldqPb); plot(1:12,qdlPb); plot(1:12,dqdlPb);
% set(gca,'XTick',1:12); xlim([1 12])
% set(gca,'XTickLabels',{}); xlabel('(e)')

% subplot(5,1,3) %export
% plot(1:12,exportPa0,1:12,exportPb0); hold all
% plot(1:12,exportPa1,'--','Color',[0 0.45 0.74])
% plot(1:12,exportPb1,'--','Color',[0.85 0.33 0.1])
% set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
% set(gca,'XTickLabel',{});%'J','F','M','A','M','J','J','A','S','O','N','DJ'}) 
% ylabel('gC/m^2 y'); xlabel('(c)'); legend('slow 2000','fast 2000','slow 2100','fast 2100')

subplot(4,1,3) %QL decomp
plot(1:12,dqlPa/QLmaxPAP(1),1:12,dqlPb/QLmaxPAP(2)); hold all;
plot(1:12,qdlPa/QLmaxPAP(1),'-o','Color',[0 0.45 0.74])
plot(1:12,ldqPb/QLmaxPAP(2),'-o','Color',[0.85 0.33 0.1])
legend('slow \Delta QL','fast \Delta QL','slow Q\Delta L','fast L\Delta Q')
set(gca,'XTick',1:12); xlim([1 12]); set(gca,'fontsize',12)
set(gca,'XTickLabels',{})
ylabel('nondim'); xlabel('(d)')

subplot(4,1,4) %xmxl
hold on
f1=fill([1:12 12:-1:1],[xmxlP([12 1:11],1).' xmxlP([11:-1:1 12],3).'],'g','LineStyle','none');
alpha(f1,0.5)
f2=fill([1:12 12:-1:1],[xmxlP1(:,1).' xmxlP1(12:-1:1,3).'],[0.5 0.5 0.5],'LineStyle','none');
alpha(f2,0.5)
plot(1:12,xmxlP(:,2),'Color',[0.1 0.6 0.1],'LineWidth',2)
plot(1:12,xmxlP1(:,2),'--','Color',[0.5 0.5 0.5],'LineWidth',2)
xlim([1 12]); set(gca,'XTick',1:12); set(gca,'fontsize',12)
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'}) 
xlabel('(e)'); ylabel('m'); legend('2000','2100')

