% Basin and Biome analyses of 12 parameter cases, now by order parameter
%% Load 12 cases production
load('globalLatlonbasin.mat')
load('yr10JN100m.mat', 'JNall2000','JNall2100','alpha','kn','mu')
taream=ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','TAREA')./1e4;%square cm to square m

prod=squeeze(mean(nansum(-JNall2000*(10*86400*365*1e-3*14*117/16).*repmat(taream,[1 1 10 12 12]),[3 1 2]),4));
prod1=squeeze(mean(nansum(-JNall2100*(10*86400*365*1e-3*14*117/16).*repmat(taream,[1 1 10 12 12]),[3 1 2]),4));

orderparam12=kn(:)./mu(:)+1./(alpha(:).*mu(:));
[orderp2,iorder]=sort(orderparam12);
orderp3=orderp2; orderp3(9)=orderp3(9)+1;  orderp3(11)=orderp3(11)+1;

%% Basin calculatios 

mask2(:,:,1)=(basin==1);%Southern
mask2(:,:,5)=(basin==2)&(lat>0);%N Pac 1
mask2(:,:,2)=(basin==2)&(lat<0);%S Pac
mask2(:,:,4)=(basin==3);%Indian
mask2(:,:,6)=(basin==6 |basin==8 |basin==9)&(lat>0)&(lat<66.5);%N Atl 
mask2(:,:,3)=(basin==6)&(lat<0);%S Atl
mask2(:,:,7)=(basin>0)&(lat>66.5);%Arct, basin 10

basin2=nan(size(lat));
for i=1:7
   basin2(mask2(:,:,i))=i; 
end

basinprod=nan(12,7);%12 cases, 7 basins
basinprod1=nan(12,7);
for i=1:7
    basinprod(:,i)=squeeze(mean(nansum(-JNall2000*(10*86400*365*1e-3*14*117/16).*repmat(taream.*mask2(:,:,i),[1 1 10 12]),[3 1 2]),4));
    basinprod1(:,i)=squeeze(mean(nansum(-JNall2100*(10*86400*365*1e-3*14*117/16).*repmat(taream.*mask2(:,:,i),[1 1 10 12]),[3 1 2]),4));
end
%% Plot Basin map, production 2000s, 2100s, % change
y1=sind(lat); pbasin=100*(basinprod1-basinprod)./basinprod;
orderp2b=round(orderp2,3,'significant');
c1=cmocean('thermo',8); colors1=zeros(7,3); colors1(5:7,:)=c1(5:7,:);
c1=cmocean('haline',8); colors1(2:4,:)=c1([2 5 7],:);
figure; subplot(2,2,1)
scatter(lon(:),y1(:),16,basin2(:),'filled'); hold on; 
scatter(lon(basin==0),y1(basin==0),9,[0.7 0.7 0.7],'filled')
xlim([0 360]); set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
cb1=colorbar; caxis([0.5 7.5]); colormap(colors1);
cb1.Ticks=1:7; 
cb1.TickLabels={'Southern','S Pacific','S Atlantic','Indian','N Pacific','N Atlantic','Arctic'};
set(gca,'fontsize',12)
title('Basins'); ylabel('latitude'); xlabel({'longitude','(a)'})
subplot(2,2,2); 
for i=1:7
plot(1:12,basinprod(iorder,i),'-o','Color',colors1(i,:)); hold on
end
xlim([0.5 12.5]); set(gca,'XTick',1:12)
set(gca,'fontsize',12)
set(gca,'XTickLabels',orderp2b)
ylabel({'new production 2000s','gC'}); xlabel({'order parameter','(b)'})
subplot(2,2,3);
for i=1:7
plot(1:12,basinprod1(iorder,i),'-o','Color',colors1(i,:)); hold on
end
xlim([0.5 12.5]); set(gca,'XTick',1:12)
set(gca,'fontsize',12)
set(gca,'XTickLabels',orderp2b)
ylabel({'new production 2100s','gC'}); xlabel({'order parameter','(c)'})
subplot(2,2,4);
for i=1:7
plot(1:12,pbasin(iorder,i),'-o','Color',colors1(i,:)); hold on
end
xlim([0.5 12.5]); set(gca,'XTick',1:12)
set(gca,'fontsize',12)
set(gca,'XTickLabels',orderp2b)
ylabel({'% change of','new production'}); xlabel({'order parameter','(d)'})
%% Biome calculations

w0=ncread('w100m2000annual.nc','WVEL');
w1=ncread('w100m2100annual.nc','WVEL');
load('stratificationNutrientQ.mat', 'lat')
load('stratificationNutrientQ.mat', 'lon')
load('stratificationNutrientQ.mat', 'xmxlmax')%cm
load('stratificationNutrientQ.mat', 'xmxlmax1')%cm
load('stratificationNutrientQ.mat', 'taream')%m^2
fice1=ncread('g.e21.G1850ECO.T62_g17.param2100.213.pop.h.IFRAC.002101-003012.nc','IFRAC');
maxfice1=nanmax(fice1(:,:,61:end),[],3);
clear fice1;
fice=ncread('g.e21.G.T62_g17.param2000.213.pop.h.IFRAC.002101-003012.nc','IFRAC');
maxfice=nanmax(fice(:,:,61:end),[],3);
clear fice;

maskc0(:,:,1)=(lat>-5)&(lat<5)&(w0<0);%equatorial downwelling
maskc0(:,:,2)=(lat>-5)&(lat<5)&(w0>0);%equatorial upwelling
maskc0(:,:,3)=((lat>5)|(lat<-5))&(xmxlmax<15000)&(w0<0);%subtrop perm strat
maskc0(:,:,4)=((lat>5)|(lat<-5))&(xmxlmax>15000)&(w0<0);%subtrop seasonal strat
maskc0(:,:,5)=(((lat>5)&(lat<30))|((lat<-5)&(lat>-35)))&(w0>0);%low lat upwelling
maskc0(:,:,6)=(maxfice1<0.1)&((lat>30)|(lat<-35))&(w0>0);%subpolar
maskc0(:,:,7)=(maxfice>0.1)&(lat>0); %N ice
maskc0(:,:,8)=(maxfice>0.1)&(lat<0); %S ice


maskc1(:,:,1)=(lat>-5)&(lat<5)&(w1<0);%equatorial downwelling
maskc1(:,:,2)=(lat>-5)&(lat<5)&(w1>0);%equatorial upwelling
maskc1(:,:,3)=((lat>5)|(lat<-5))&(xmxlmax1<15000)&(w1<0);%subtrop perm strat
maskc1(:,:,4)=((lat>5)|(lat<-5))&(xmxlmax1>15000)&(w1<0);%subtrop seasonal strat
maskc1(:,:,5)=(((lat>5)&(lat<30))|((lat<-5)&(lat>-35)))&(w1>0);%low lat upwelling
maskc1(:,:,6)=(maxfice1<0.1)&((lat>30)|(lat<-35))&(w1>0);%subpolar
maskc1(:,:,7)=(maxfice1>0.1)&(lat>0); %N ice
maskc1(:,:,8)=(maxfice1>0.1)&(lat<0); %S ice

biomecarea=squeeze(nansum(nansum(maskc0.*repmat(taream,[1 1 8]))));
biomecarea1=squeeze(nansum(nansum(maskc1.*repmat(taream,[1 1 8]))));

biomec=nan(size(lat));
biomec1=nan(size(lat));
for i=1:8
   biomec(maskc0(:,:,i))=i; 
   biomec1(maskc1(:,:,i))=i; 
end

biomecprod=nan(12,8);%12 cases, 12 biomes
biomecprod1=nan(12,8);
for i=1:8
    biomecprod(:,i)=squeeze(mean(nansum(-JNall2000*(10*86400*365*1e-3*14*117/16).*repmat(taream.*maskc0(:,:,i),[1 1 10 12 12]),[3 2 1]),4));
    biomecprod1(:,i)=squeeze(mean(nansum(-JNall2100*(10*86400*365*1e-3*14*117/16).*repmat(taream.*maskc1(:,:,i),[1 1 10 12 12]),[3 2 1]),4));
end

%% Biome plot: 
%map2000, map2100
%areas2000,%change
%prodrate2000, %change
%totprod2000, %change

y1=sind(lat);  colors1=zeros(8,3);
c1=cmocean('thermo',8); colors1(5:7,:)=c1(5:7,:);
c2=cmocean('haline',8); colors1([1 3 4],:)=c2([2 5 7],:); colors1(2,:)=[0.9 0 0.1];

figure; subplot(4,2,1)
scatter(lon(:),y1(:),16,biomec(:),'filled'); hold on; scatter(lon(basin==0),y1(basin==0),9,[0.7 0.7 0.7],'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'fontsize',12)
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
ylabel('latitude'); xlabel({'longitude','(a)'})
xlim([-180 180])
subplot(4,2,2)
scatter(lon(:),y1(:),16,biomec1(:),'filled'); hold on; scatter(lon(basin==0),y1(basin==0),9,[0.7 0.7 0.7],'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'fontsize',12)
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
cb2=colorbar; cb2.Ticks=1:8; caxis([0.5 8.5]);
cb2.TickLabels={'Eq D','Eq U','ST PS','ST SS','LLU','SP','N I','S I'};
 xlabel({'longitude','(b)'}); colormap(colors1); xlim([-180 180])
 
 subplot(4,2,3)
for i=1:8; hold on; plot(i,biomecarea(i),'o','MarkerSize',10,'Color',colors1(i,:),'MarkerFaceColor',colors1(i,:)); end
%for i=1:8; hold on; plot(i,biomecarea1(i),'d','MarkerSize',10,'Color',colors1(i,:),'MarkerFaceColor',colors1(i,:)); end
set(gca,'fontsize',12)
xlim([0.5 8.5]); set(gca,'XTick',1:8);
set(gca,'XTickLabels',{'Eq D','Eq U','ST PS','ST SS','LLU','SP','N I','S I'});
ylabel('biome area, m^2'); xlabel({'biome','(c)'})

subplot(4,2,4)
parea=100*(biomecarea1-biomecarea)./biomecarea;
for i=1:8; hold on; plot(i,parea(i),'s','MarkerSize',10,'Color',colors1(i,:),'MarkerFaceColor',colors1(i,:)); end
plot([0 9],[0 0],'k')
set(gca,'fontsize',12)
xlim([0.5 8.5]); set(gca,'XTick',1:8);
set(gca,'XTickLabels',{'Eq D','Eq U','ST PS','ST SS','LLU','SP','N I','S I'});
ylabel({'biome area','% change'}); xlabel({'biome','(d)'})

subplot(4,2,5)
for i=1:8; hold on; plot(1:12,biomecprod(:,i)./biomecarea(i),'-o','Color',colors1(i,:)); end
set(gca,'fontsize',12)
xlim([0.5 12.5]); set(gca,'XTick',1:12); set(gca,'XTickLabels',orderp2b)
xlabel({'order parameter','(e)'}); ylabel({'production rate','gC/m^2'})

subplot(4,2,6) 
prate=100*((biomecprod1./repmat(biomecarea1.',[12 1]))-(biomecprod./repmat(biomecarea.',[12 1])))./(biomecprod./repmat(biomecarea.',[12 1]));
for i=1:8; hold on; plot(1:12,prate(:,i),'-o','Color',colors1(i,:)); end
set(gca,'fontsize',12)
xlim([0.5 12.5]); set(gca,'XTick',1:12); set(gca,'XTickLabels',orderp2b)
xlabel({'order parameter','(e)'}); ylabel({'production rate','% change'})

subplot(4,2,7)
for i=1:8; hold on; plot(1:12,biomecprod(:,i),'-o','Color',colors1(i,:)); end
set(gca,'fontsize',12)
xlim([0.5 12.5]); set(gca,'XTick',1:12); set(gca,'XTickLabels',orderp2b)
xlabel({'order parameter','(e)'}); ylabel({'annual production','gC'})

subplot(4,2,8) 
prate2=100*(biomecprod1-biomecprod)./(biomecprod);
for i=1:8; hold on; plot(1:12,prate2(:,i),'-o','Color',colors1(i,:)); end
set(gca,'fontsize',12)
xlim([0.5 12.5]); set(gca,'XTick',1:12); set(gca,'XTickLabels',orderp2b)
xlabel({'order parameter','(e)'}); ylabel({'annual production','% change'})
%% Plot Biome area, frac change, production rate, frac change, tot prod, frac change
%old versions
figure; subplot(4,1,1)
hold all; plot(1:8,biomecprod([1:3 5:10 12],:).'./repmat(biomecarea,[1 10]),'ko','MarkerSize',4); plot(1:8,biomecprod([4 11],:).'./repmat(biomecarea,[1 2]),'*','MarkerSize',8);
xlim([0.5 8.5]); set(gca,'XTick',1:8); set(gca,'XTickLabels',{}); set(gca,'fontsize',12)
ylabel({'Production 2000s','gC/m^2 yr'}); ylim([0 90])
subplot(4,1,2)
hold all; plot(1:8,biomecprod1(4,:).'./biomecarea1,'*','MarkerSize',8,'Color',[0.3 0.75 0.93]); plot(1:8,biomecprod1(11,:).'./biomecarea1,'*','MarkerSize',8);
plot(1:8,biomecprod1([1:3 5:10 12],:).'./repmat(biomecarea1,[1 10]),'ko','MarkerSize',4); 
plot(1:8,biomecprod1(4,:).'./biomecarea1,'*'); plot(1:8,biomecprod1(11,:).'./biomecarea1,'*','MarkerSize',8);
legend('slow','fast','10 other cases'); ylim([0 90])
xlim([0.5 8.5]); set(gca,'XTick',1:8); set(gca,'XTickLabels',{}); set(gca,'fontsize',12)
ylabel({'Production 2100s','gC/m^2 yr'})
subplot(4,1,3)
hold all; plot(1:8,biomecprod1([1:3 5:10 12],:).'./biomecprod([1:3 5:10 12],:).','ko','MarkerSize',4);plot(1:8,biomecprod1([4 11],:).'./biomecprod([4 11],:).','*','MarkerSize',8); 
xlim([0.5 8.5]); set(gca,'XTick',1:8); set(gca,'XTickLabels',{}); set(gca,'fontsize',12)
ylabel({'Production 2100/2000','nondim'})
subplot(4,1,4)
hold all; plot(1:8,biomecarea1./biomecarea,'rd'); plot([0 9],[1 1],'k')
xlim([0.5 8.5]); set(gca,'XTick',1:8);  set(gca,'fontsize',12)
ylabel({'Biome area 2100/2000','nondim'}); set(gca,'XTickLabels',{'Eq D','Eq U','ST PS','ST SS','LLU','SP','N I','S I'})


figure; subplot(4,1,1)
hold all; plot(1:8,biomecprod([1:3 5:10 12],:).'./repmat(biomecarea,[1 10]),'ko','MarkerSize',4); plot(1:8,biomecprod([4 11],:).'./repmat(biomecarea,[1 2]),'*','MarkerSize',8);
xlim([0.5 8.5]); set(gca,'XTick',1:8); set(gca,'XTickLabels',{}); set(gca,'fontsize',12)
ylabel({'Production rate 2000s','gC/m^2 yr'}); ylim([0 90])
subplot(4,1,2)
hold all; plot(1:8,biomecprod1(4,:).'.*biomecarea./(biomecprod(4,:).'.*biomecarea1),'*','MarkerSize',8,'Color',[0.3 0.75 0.93]); plot(1:8,biomecprod1(11,:).'.*biomecarea./(biomecprod(11,:).'.*biomecarea1),'*','MarkerSize',8);
plot(1:8,biomecprod1([1:3 5:10 12],:).'.*repmat(biomecarea,[1 10])./(biomecprod([1:3 5:10 12],:).'.*repmat(biomecarea1,[1 10])),'ko','MarkerSize',4); 
plot(1:8,biomecprod1(4,:).'.*biomecarea./(biomecprod(4,:).'.*biomecarea1),'*'); plot(1:8,biomecprod1(11,:).'.*biomecarea./(biomecprod(11,:).'.*biomecarea1),'*','MarkerSize',8);
plot([0 9],[1 1],'k')
legend('slow','fast','10 other cases'); ylim([0 90])
xlim([0.5 8.5]); set(gca,'XTick',1:8); set(gca,'XTickLabels',{}); set(gca,'fontsize',12)
ylabel({'Production rate 2100/2000','nondim'})
subplot(4,1,3)
hold all; plot(1:8,biomecprod1([1:3 5:10 12],:).'./biomecprod([1:3 5:10 12],:).','ko','MarkerSize',4);plot(1:8,biomecprod1([4 11],:).'./biomecprod([4 11],:).','*','MarkerSize',8); 
xlim([0.5 8.5]); set(gca,'XTick',1:8); set(gca,'XTickLabels',{}); set(gca,'fontsize',12)
ylabel({'Production 2100/2000','gC/gC'})
subplot(4,1,4)
hold all; plot(1:8,biomecarea1./biomecarea,'rd'); plot([0 9],[1 1],'k')
xlim([0.5 8.5]); set(gca,'XTick',1:8);  set(gca,'fontsize',12)
ylabel({'Biome area 2100/2000','m^2/m^2'}); set(gca,'XTickLabels',{'Eq D','Eq U','ST PS','ST SS','LLU','SP','N I','S I'})

figure;subplot(1,2,1)
bar(biomecprod,'stacked')
xlabel('various parameters')
title('Biome production 2000')
ylabel('New production, gC/yr')
subplot(1,2,2)
bar(biomecprod1,'stacked')
title('Biome production 2000')
legend('Eq D','Eq U','ST PS','ST SS','LLU','SP','N I','S I')


