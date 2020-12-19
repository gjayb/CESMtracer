
%% load N, 100m-average
load('yr10N100m.mat', 'Nall2000','Nall2100')
N=squeeze(mean(Nall2000(:,:,:,:,[4 11]),3));
N1=squeeze(mean(Nall2100(:,:,:,:,[4 11]),3));
clear Nall2*
taream=ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','TAREA')./1e4;
load('climateAndRegions.mat', 'subtropSPac','osmosis')%,'arctic')
load('globalLatlonbasin.mat')
arctic=(lat>66.5)&(basin~=0);
%% seasonal cycles: globe, ssp, arctic, pap
logN=(lat>0)&(basin~=0);%northern and southern hemisphere
logS=(lat<0)&(basin~=0);
arean=sum(taream(logN));
areas=sum(taream(logS));
areat=sum(taream(basin~=0));
areaA=sum(taream(arctic));
areaS=sum(taream(subtropSPac));
areaP=sum(taream(osmosis));

for j=1:2
    for i=1:12
        hold1=N(:,:,i,j).*taream;
        Nn(i,j)=nansum(hold1(logN));
        hold1=N1(:,:,i,j).*taream;
        Nn1(i,j)=nansum(hold1(logN));
        hold1=N(:,:,i,j).*taream;
        Ns(i,j)=nansum(hold1(logS));
        hold1=N1(:,:,i,j).*taream;
        Ns1(i,j)=nansum(hold1(logS));
        hold1=N(:,:,i,j).*taream;
        Na(i,j)=nansum(hold1(arctic))/areaA;
        hold1=N1(:,:,i,j).*taream;
        Na1(i,j)=nansum(hold1(arctic))/areaA;
        hold1=N(:,:,i,j).*taream;
        Nsp(i,j)=nansum(hold1(subtropSPac))/areaS;
        hold1=N1(:,:,i,j).*taream;
        Nsp1(i,j)=nansum(hold1(subtropSPac))/areaS;
        hold1=N(:,:,i,j).*taream;
        Np(i,j)=nansum(hold1(osmosis))/areaP;
        hold1=N1(:,:,i,j).*taream;
        Np1(i,j)=nansum(hold1(osmosis))/areaP;
    end
end

Ng=(Nn+Ns([7:12 1:6],:))/areat; Ng1=(Nn1+Ns1([7:12 1:6],:))/areat;
%% load WOA seasonal cycles for all regions
%pap 11-27W, 40-52N
%SSP 10-35S, 143-287E (73-216W)
%arcic >66.5N
load('woaNitrateSeasonal2.mat')
%% plot seasonal cycles

figure; 
subplot(2,2,1) %global
%f1=fill([1:12 12:-1:1],[nitrate_ten_glob(:,1).' nitrate_ten_glob(12:-1:1,2).'],[0.7 0.7 0.7],'LineStyle','none');
%alpha(f1,0.4)
hold on; plot(1:12,nitrate_glob,'k--','LineWidth',2)
plot(1:12,Ng([12 1:11],:),1:12,Ng1([12 1:11],:),'--')
set(gca,'XTick',1:12); xlim([1 12]);
set(gca,'XTickLabels',{}); ylabel('mmolN/m^3'); set(gca,'fontsize',12)
xlabel('(a)')

subplot(2,2,2)%ssp
%f1=fill([1:12 12:-1:1],[nitrate_ten_ssp(:,1).' nitrate_ten_ssp(12:-1:1,2).'],[0.7 0.7 0.7],'LineStyle','none');
%alpha(f1,0.4)
hold on; plot(1:12,nitrate_ssp,'k--','LineWidth',2)
plot(1:12,Nsp([12 1:11],:),1:12,Nsp1([12 1:11],:),'--')
set(gca,'XTick',1:12); xlim([1 12]);
set(gca,'XTickLabels',{}); set(gca,'fontsize',12)
xlabel('(b)')

subplot(2,2,3)%arctic
%f1=fill([1:12 12:-1:1],[nitrate_ten_arctic(:,1).' nitrate_ten_arctic(12:-1:1,2).'],[0.7 0.7 0.7],'LineStyle','none');
%alpha(f1,0.4)
hold on; plot(1:12,nitrate_arctic,'k--','LineWidth',2)
plot(1:12,Na([12 1:11],:),1:12,Na1([12 1:11],:),'--')
set(gca,'XTick',1:12); xlim([1 12]);
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'}) 
ylabel('mmolN/m^3'); set(gca,'fontsize',12)
xlabel('(c)')

subplot(2,2,4) %pap/osmosis
%f1=fill([1:12 12:-1:1],[nitrate_ten_osmosis(:,1).' nitrate_ten_osmosis(12:-1:1,2).'],[0.7 0.7 0.7],'LineStyle','none');
%alpha(f1,0.4)
hold on; plot(1:12,nitrate_osmosis,'k--','LineWidth',2)
plot(1:12,Np([12 1:11],:),1:12,Np1([12 1:11],:),'--')
set(gca,'XTick',1:12); xlim([1 12]);
set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'}) 
set(gca,'fontsize',12); xlabel('(d)')
legend('mean WOA nitrate','2000 slow','2000 fast','2100 slow','2100 fast')


