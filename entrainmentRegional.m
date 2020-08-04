%entrainment N fluxes breakdown
%% load budget terms
load('npFluxes100m2100tau30w10.332.mat', 'Wnutriday','DiaNutriday','DiffNutriday','taream')
%Dia from DIA_IMPVF_NUTRI, implicit diffusion, should be kpp
%Diff from HDIFB_NUTRI, should be GM-Redi
%units! these are fluxes mmol per day, should multiply by 0.365*14*117/16
%for gC/yr, and dividing by taream for gC/m^2yr
wnb1=mean(Wnutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
kppb1=mean(DiaNutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
gmrb1=mean(DiffNutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
load('npFluxes100m2000tau30w10.332.mat', 'Wnutriday','DiaNutriday','DiffNutriday')
wnb=mean(Wnutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
kppb=mean(DiaNutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
gmrb=mean(DiffNutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
load('npFluxes100m2000tau1yw5.021.mat', 'Wnutriday','DiaNutriday','DiffNutriday')
wna=mean(Wnutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
kppa=mean(DiaNutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
gmra=mean(DiffNutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
load('npFluxes100m2100tau1yw5.021.mat', 'Wnutriday','DiaNutriday','DiffNutriday')
wna1=mean(Wnutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
kppa1=mean(DiaNutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;
gmra1=mean(DiffNutriday(:,:,end-11:end),3)*0.365*14*117/16./taream;

load('climateAndRegions.mat', 'subtropSPac','osmosis','arctic')

%% global annual means
globalwnb1=squeeze(areaweightedmean(areaweightedmean(wnb1,taream,2),mean(taream,2),1));
globalwnb=squeeze(areaweightedmean(areaweightedmean(wnb,taream,2),mean(taream,2),1));
globalwna1=squeeze(areaweightedmean(areaweightedmean(wna1,taream,2),mean(taream,2),1));
globalwna=squeeze(areaweightedmean(areaweightedmean(wna,taream,2),mean(taream,2),1));

globalkppb1=squeeze(areaweightedmean(areaweightedmean(kppb1,taream,2),mean(taream,2),1));
globalkppb=squeeze(areaweightedmean(areaweightedmean(kppb,taream,2),mean(taream,2),1));
globalkppa1=squeeze(areaweightedmean(areaweightedmean(kppa1,taream,2),mean(taream,2),1));
globalkppa=squeeze(areaweightedmean(areaweightedmean(kppa,taream,2),mean(taream,2),1));

globalgmrb1=squeeze(areaweightedmean(areaweightedmean(gmrb1,taream,2),mean(taream,2),1));
globalgmrb=squeeze(areaweightedmean(areaweightedmean(gmrb,taream,2),mean(taream,2),1));
globalgmra1=squeeze(areaweightedmean(areaweightedmean(gmra1,taream,2),mean(taream,2),1));
globalgmra=squeeze(areaweightedmean(areaweightedmean(gmra,taream,2),mean(taream,2),1));

figure; subplot(1,2,1)
bar([globalkppa,globalkppa1,globalkppa1-globalkppa;globalgmra,globalgmra1,globalgmra1-globalgmra;globalwna,globalwna1,globalwna1-globalwna])
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'KPP','GM-Redi','Adv'}); title('slow'); ylabel('gC/m^2 yr')
subplot(1,2,2)
bar([globalkppb,globalkppb1,globalkppb1-globalkppb;globalgmrb,globalgmrb1,globalgmrb1-globalgmrb;globalwnb,globalwnb1,globalwnb1-globalwnb])
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'KPP','GM-Redi','Adv'}); title('fast'); legend('2000','2100','\Delta')

%% SSP annual means
sspwnb1=squeeze(areaweightedmean(areaweightedmean(wnb1.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));
sspwnb=squeeze(areaweightedmean(areaweightedmean(wnb.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));
sspwna1=squeeze(areaweightedmean(areaweightedmean(wna1.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));
sspwna=squeeze(areaweightedmean(areaweightedmean(wna.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));

sspkppb1=squeeze(areaweightedmean(areaweightedmean(kppb1.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));
sspkppb=squeeze(areaweightedmean(areaweightedmean(kppb.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));
sspkppa1=squeeze(areaweightedmean(areaweightedmean(kppa1.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));
sspkppa=squeeze(areaweightedmean(areaweightedmean(kppa.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));

sspgmrb1=squeeze(areaweightedmean(areaweightedmean(gmrb1.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));
sspgmrb=squeeze(areaweightedmean(areaweightedmean(gmrb.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));
sspgmra1=squeeze(areaweightedmean(areaweightedmean(gmra1.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));
sspgmra=squeeze(areaweightedmean(areaweightedmean(gmra.*subtropSPac,taream.*subtropSPac,2),mean(taream.*subtropSPac,2),1));

figure; subplot(1,3,1)
bar([sspkppa,sspkppa1,sspkppa1-sspkppa;sspgmra,sspgmra1,sspgmra1-sspgmra;sspwna,sspwna1,sspwna1-sspwna])
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'KPP','GM-Redi','Adv'}); title('slow'); ylabel('gC/m^2 yr')
set(gca,'fontsize',12); xlabel('(a)')
subplot(1,3,2)
bar([sspkppb,sspkppb1,sspkppb1-sspkppb;sspgmrb,sspgmrb1,sspgmrb1-sspgmrb;sspwnb,sspwnb1,sspwnb1-sspwnb])
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'KPP','GM-Redi','Adv'}); title('fast'); legend('2000','2100','\Delta')
set(gca,'fontsize',12);xlabel('(b)')

% compute N profiles
load('stratificationNutrientQ.mat', 'N*')
Na0=mean(Na0,4);
Na1=mean(Na1,4);
Nb0=mean(Nb0,4);
Nb1=mean(Nb1,4);
zm=ncread('salt2000yr30.nc','z_t')/100;
%%
for i=1:60
   hold0=Na0(:,:,i);
   hold1=Na1(:,:,i);
   dNa(i)=nansum(hold1(subtropSPac).*taream(subtropSPac))/sum(taream(subtropSPac)) -nansum(hold0(subtropSPac).*taream(subtropSPac))/sum(taream(subtropSPac));
   hold0=Nb0(:,:,i);
   hold1=Nb1(:,:,i);
   dNb(i)=nansum(hold1(subtropSPac).*taream(subtropSPac))/sum(taream(subtropSPac)) -nansum(hold0(subtropSPac).*taream(subtropSPac))/sum(taream(subtropSPac));
end

subplot(1,3,3)
plot(dNa,-zm,dNb,-zm); 
set(gca,'fontsize',12);
legend('slow','fast'); xlabel({'N','(c)'}); ylabel('depth (m)')
ylim([-1000 0])
%% Arctic annual means 
arcwnb1=squeeze(areaweightedmean(areaweightedmean(wnb1.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));
arcwnb=squeeze(areaweightedmean(areaweightedmean(wnb.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));
arcwna1=squeeze(areaweightedmean(areaweightedmean(wna1.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));
arcwna=squeeze(areaweightedmean(areaweightedmean(wna.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));

arckppb1=squeeze(areaweightedmean(areaweightedmean(kppb1.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));
arckppb=squeeze(areaweightedmean(areaweightedmean(kppb.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));
arckppa1=squeeze(areaweightedmean(areaweightedmean(kppa1.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));
arckppa=squeeze(areaweightedmean(areaweightedmean(kppa.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));

arcgmrb1=squeeze(areaweightedmean(areaweightedmean(gmrb1.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));
arcgmrb=squeeze(areaweightedmean(areaweightedmean(gmrb.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));
arcgmra1=squeeze(areaweightedmean(areaweightedmean(gmra1.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));
arcgmra=squeeze(areaweightedmean(areaweightedmean(gmra.*arctic,taream.*arctic,2),mean(taream.*arctic,2),1));

figure; subplot(1,2,1)
bar([arckppa,arckppa1,arckppa1-arckppa;arcgmra,arcgmra1,arcgmra1-arcgmra;arcwna,arcwna1,arcwna1-arcwna])
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'KPP','GM-Redi','Adv'}); title('slow'); ylabel('gC/m^2 yr')
subplot(1,2,2)
bar([arckppb,arckppb1,arckppb1-arckppb;arcgmrb,arcgmrb1,arcgmrb1-arcgmrb;arcwnb,arcwnb1,arcwnb1-arcwnb])
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'KPP','GM-Redi','Adv'}); title('fast'); legend('2000','2100','\Delta')

%% PAP annual means 

papwnb1=squeeze(areaweightedmean(areaweightedmean(wnb1.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));
papwnb=squeeze(areaweightedmean(areaweightedmean(wnb.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));
papwna1=squeeze(areaweightedmean(areaweightedmean(wna1.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));
papwna=squeeze(areaweightedmean(areaweightedmean(wna.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));

papkppb1=squeeze(areaweightedmean(areaweightedmean(kppb1.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));
papkppb=squeeze(areaweightedmean(areaweightedmean(kppb.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));
papkppa1=squeeze(areaweightedmean(areaweightedmean(kppa1.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));
papkppa=squeeze(areaweightedmean(areaweightedmean(kppa.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));

papgmrb1=squeeze(areaweightedmean(areaweightedmean(gmrb1.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));
papgmrb=squeeze(areaweightedmean(areaweightedmean(gmrb.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));
papgmra1=squeeze(areaweightedmean(areaweightedmean(gmra1.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));
papgmra=squeeze(areaweightedmean(areaweightedmean(gmra.*osmosis,taream.*osmosis,2),mean(taream.*osmosis,2),1));

figure; subplot(1,2,1)
bar([papkppa,papkppa1,papkppa1-papkppa;papgmra,papgmra1,papgmra1-papgmra;papwna,papwna1,papwna1-papwna])
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'KPP','GM-Redi','Adv'}); title('slow'); ylabel('gC/m^2 yr')
subplot(1,2,2)
bar([papkppb,papkppb1,papkppb1-papkppb;papgmrb,papgmrb1,papgmrb1-papgmrb;papwnb,papwnb1,papwnb1-papwnb])
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'KPP','GM-Redi','Adv'}); title('fast'); legend('2000','2100','\Delta')


