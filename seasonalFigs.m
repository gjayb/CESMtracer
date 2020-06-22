% Seasonal cycle plots for globe and 3 regions

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

%load XMXL %probably not using this?
XMXL=0.01*(ncread('g.e21.G.T62_g17.param2000.121.pop.h.XMXL.002101-003012.nc','XMXL',[1 1 109],[Inf Inf Inf]));%0.01 to go cm-->m
XMXL1=0.01*(ncread('g.e21.G1850ECO.T62_g17.param2100.121.pop.h.XMXL.002101-003012.nc','XMXL',[1 1 109],[Inf Inf Inf]));
dXMXL=XMXL1-XMXL;

%load stratification
PD=1000*(ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','PD',[1 1 1 109],[Inf Inf 16 Inf])-1);
PD1=1000*(ncread('g.e21.G1850ECO.T62_g17.param2100.123.pop.h.PD.002101-003012.nc','PD',[1 1 1 109],[Inf Inf 16 Inf])-1);
z=ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','z_t')./100;
strat=squeeze(PD(:,:,16,end-11:end)-PD(:,:,6,end-11:end))./100;%year 10  potential density at 155m-55m/100m
strat1=squeeze(PD1(:,:,16,end-11:end)-PD1(:,:,6,end-11:end))./100;
dstrat=strat1-strat;
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
clear Q0 Q1 L0 L1

%% global 6-month offset cycle computation

logN=(lat>0)&(basin~=0);%northern and southern hemisphere
logS=(lat<0)&(basin~=0);
arean=sum(taream(logN));
areas=sum(taream(logS));
areat=sum(taream(basin~=0));

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
entrainGb0=((entrainBn([36 25:35])+entrainBn([24 13:23])+entrainBn([12 1:11]))+(entrainBs([30:36 25:29])+entrainBs([24 13:23])+entrainBs([12 1:11])))*0.365*14*117/(16*3*areat);
entrainGb1=((entrainBn1([36 25:35])+entrainBn1([24 13:23])+entrainBn1([12 1:11]))+(entrainBs1([30:36 25:29])+entrainBs1([24 13:23])+entrainBs1([12 1:11])))*0.365*14*117/(16*3*areat);
entrainGa0=((entrainAn([36 25:35])+entrainAn([24 13:23])+entrainAn([12 1:11]))+(entrainAs([30:36 25:29])+entrainAs([24 13:23])+entrainAs([12 1:11])))*0.365*14*117/(16*3*areat);
entrainGa1=((entrainAn1([36 25:35])+entrainAn1([24 13:23])+entrainAn1([12 1:11]))+(entrainAs1([30:36 25:29])+entrainAs1([24 13:23])+entrainAs1([12 1:11])))*0.365*14*117/(16*3*areat);

stratG1=((stratn1([12 1:11]))+(strats1([12 1:11])))/(3*areat);
stratG0=((stratn([12 1:11]))+(strats([12 1:11])))/(3*areat);
mldG1=((mldn1([36 25:35])+mldn1([24 13:23])+mldn1([12 1:11]))+(mlds1([30:36 25:29])+mlds1([24 13:23])+mlds1([12 1:11])))/(3*areat);
mldG0=((mldn([36 25:35])+mldn([24 13:23])+mldn([12 1:11]))+(mlds([30:36 25:29])+mlds([24 13:23])+mlds([12 1:11])))/(3*areat);

QGa0=(QAn([12 1:11])+QAs([12 1:11]))/areat;
QGb0=(QBn([12 1:11])+QBs([12 1:11]))/areat;
QGa1=(QAn1([12 1:11])+QAs1([12 1:11]))/areat;
QGb1=(QBn1([12 1:11])+QBs1([12 1:11]))/areat;

LGa0=(LAn([12 1:11])+LAs([12 1:11]))/areat;
LGb0=(LBn([12 1:11])+LBs([12 1:11]))/areat;
LGa1=(LAn1([12 1:11])+LAs1([12 1:11]))/areat;
LGb1=(LBn1([12 1:11])+LBs1([12 1:11]))/areat;
%% global plot
%entrain,prod,export; strat,mld; Q1/Q,L1/L

figure;
subplot(3,2,1) %entrain

subplot(3,2,2) %Q1/Q, L1/L

subplot(3,2,3) %prod

subplot(3,2,4) %mld

subplot(3,2,5) %export

subplot(3,2,6) %strat


%% subtropical south pacific computation
%logSSP=;
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

%% SSP plot
% Nflux, production, P flux down, growth limit changes
figure;
subplot(4,1,1) %entrain

subplot(4,1,2) %prod

subplot(4,1,3) %export

subplot(4,1,4) %Q1/Q, L1/L


%% Arctic computation
% Nflux, production, P flux down, growth limits, incoming light

%% arctic plot

figure;
subplot(3,2,1) %entrain

subplot(3,2,2) %Q

subplot(3,2,3) %prod

subplot(3,2,4) %L

subplot(3,2,5) %export

subplot(3,2,6) %light

%% Porcupine Abyssal Plane computation
% Nflux, production, P flux down, growth limit changes, MLD as XMXL
maskosmosis=(lon>-26)&(lon<-16)&(lat>40)&(lat<50);

%% PAP plot

figure;
subplot(3,2,1) %entrain

subplot(3,2,2) %Q1/Q, L1/L

subplot(3,2,3) %prod

subplot(3,2,4) %XMXL

subplot(3,2,5) %export

