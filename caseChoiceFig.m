%Figure 2:
%annual zonal surface mean plots of N, P for 27 cases, and surface N vs
%case timescale (order parameter)
%Figure 3: annual total new production, surface N in both climates for 12 cases

%load zonal annual means
load('choosingParamsNnew.mat', 'zonalNO3woa','zonalN','zonalNb','zonalTarea','latMean','alphas','mus','kns')
load('choosingParamsPnew.mat','zonalP','zonalBEC_C')
zonalPb=reshape(zonalP,[58 60 27]);
better12=zeros([3 3 3]);
better12(1,2,:)=1; better12(2,1,:)=1; better12(3,1,:)=1; better12(3,3,:)=1;
%better12 [2:4 9 11:13 18 20:22 27]
slowcase=zeros([3 3 3]); slowcase(1,2,1)=1; %case 4
fastcase=zeros([3 3 3]); fastcase(3,3,2)=1; %case 18

alpha=repmat(alphas,[1 3 3]);
kn(1,1,1:3)=kns; kn=repmat(kn,[3 3 1]);
mu(1,1:3,1)=mus; mu=repmat(mu,[3 1 3]);
%% plot zonal annual means
dark1=logical(reshape(better12,27,1));
figure; subplot(1,3,1); 
plot(zonalNO3woa(:,1),latMean,'k','Linewidth',2); hold on
plot(zonalNb(:,1,4),latMean,'Color',[0 0.45 0.74],'Linewidth',2); plot(zonalNb(:,1,18),latMean,'Color',[0.85 0.33 0.1],'Linewidth',2);
plot(squeeze(zonalNb(:,1,[1:3 5:17 19:27])),latMean,'Color',[0.6 0.6 0.6]); 
plot(zonalNO3woa(:,1),latMean,'k','Linewidth',2); hold on
plot(zonalNb(:,1,4),latMean,'Color',[0 0.45 0.74],'Linewidth',2); plot(zonalNb(:,1,18),latMean,'Color',[0.85 0.33 0.1],'Linewidth',2);
legend('WOA NO_3','slow','fast','25 other cases'); ylim([-80 90]); set(gca,'fontsize',12)
ylabel('latitude'); xlabel({'N units','(a)'})

subplot(1,3,3) %redfield in BEC is 16N:117C
plot(zonalBEC_C(:,1).*16./117,latMean,'k','Linewidth',2); hold on
plot(zonalPb(:,1,4),latMean,'Color',[0 0.45 0.74],'Linewidth',2); plot(zonalPb(:,1,18),latMean,'Color',[0.85 0.33 0.1],'Linewidth',2);
plot(squeeze(zonalPb(:,1,2)),latMean,'Color',[0.3 0.3 0.3]); 
plot(squeeze(zonalPb(:,1,1)),latMean,'Color',[0.8 0.8 0.8]); 
plot(squeeze(zonalPb(:,1,[5:8 10 14:17 19 23:26])),latMean,'Color',[0.8 0.8 0.8]); 
plot(squeeze(zonalPb(:,1,[3:4 9 11:13 18 20:22 27])),latMean,'Color',[0.3 0.3 0.3]); 
plot(zonalBEC_C(:,1).*16./117,latMean,'k','Linewidth',2); hold on
plot(zonalPb(:,1,4),latMean,'Linewidth',2,'Color',[0 0.45 0.74]); plot(zonalPb(:,1,18),latMean,'Linewidth',2,'Color',[0.85 0.33 0.1]);
legend('BEC','slow','fast','10 next best cases','15 other cases'); ylim([-80 90]); set(gca,'fontsize',12)
xlabel({'N units','(b)'}); ylabel('latitude');%set(gca,'YTickLabels',{})
%% add param summary value/ordering at max N
orderparam=kn(:)./mu(:)+1./(alpha(:).*mu(:));%kn(:)./mu(:)+1./(alpha(:).*mu(:));
for i=1:27
   %[xval(i),lati]=max(zonalNb(:,1,i)); 
   %yval(i)=latMean(lati);
   j=mod(i,6)+1;
   xval(i)=zonalNb(j,1,i);
   yval(i)=latMean(j);
end
%subplot(2,3,1); hold on
%scatter(xval,yval,16,orderparam,'filled'); colorbar
%set(gca,'ColorScale','log'); legend('WOA NO_3','slow','fast','25 other cases','order parameter');
subplot(1,3,2)
Nsurf=nansum(squeeze(zonalNb(:,1,:)).*repmat(zonalTarea(:,1),[1 27]),1)./nansum(repmat(zonalTarea(:,1),[1 27]),1);
plot(Nsurf,orderparam,'ko')
hold on; plot(Nsurf(4),orderparam(4),'*','Color',[0 0.45 0.74])
plot(Nsurf(18),orderparam(18),'*','Color',[0.85 0.33 0.1])
set(gca,'fontsize',12); xlabel({'N','(b)'}); ylabel('case timescale, days')
%% load 12 cases production
load('yr10JN100m.mat', 'JNall2000','JNall2100','alpha','kn','mu')
taream=ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','TAREA')./1e4;%square cm to square m
%Jnutri100day=squeeze(sum(repmat(dz2,[320 384 1 120]).*repmat(taream,[1 1 10 120]).*Jnutri(:,:,1:10,:)*86400*365*1e-3*14,3));
%J_NUTRI is the rate mmol/(s m^3) of transfer N-->P
%convert to PgC/yr: J_NUTRI*(dz*taream m^3)*(86400*365 s/yr)*(1e-3mol/mmol)*(14gN/mol)*(117C/16N)
prod=squeeze(mean(nansum(-JNall2000*(10*86400*365*1e-3*14*117/16).*repmat(taream,[1 1 10 12 12]),[3 1 2]),4));
prod1=squeeze(mean(nansum(-JNall2100*(10*86400*365*1e-3*14*117/16).*repmat(taream,[1 1 10 12 12]),[3 1 2]),4));
clear JNall2*

orderparam12=kn(:)./mu(:)+1./(alpha(:).*mu(:));
[orderp2,iorder]=sort(orderparam12);
orderp3=orderp2; orderp3(9)=orderp3(9)+1;  orderp3(11)=orderp3(11)+1;
%% plot 12 cases production and % change
figure; subplot(2,2,1)
plot(1:12,prod(iorder),'ko'); hold on; plot(1:12,prod1(iorder),'rd')
set(gca,'XTick',1:12); set(gca,'XTickLabels',{})
ylabel('new production, gC'); xlabel('(a)')
xlim([0.5 12.5]); ylim([4.5e15 8e15])
set(gca,'fontsize',12); legend('2000','2100')
subplot(2,2,3)
plot(1:12,100*(prod1(iorder)-prod(iorder))./prod(iorder),'k*')
hold on; plot(2,100*(prod1(11)-prod(11))./prod(11),'o','Color',[0.85 0.33 0.1],'MarkerSize',9)
plot(9,100*(prod1(4)-prod(4))./prod(4),'o','Color',[0 0.45 0.74],'MarkerSize',9)
ylabel('% change')
set(gca,'XTick',1:12); set(gca,'XTickLabels',num2str(round(orderp2,1)))
xlabel({'case timescale','(b)'})
xlim([0.5 12.5]); ylim([-20 0])
set(gca,'fontsize',12)

%% load 12 cases N, compute surface values
load('yr10N100m.mat', 'Nall2000','Nall2100')
Nall2000=squeeze(mean(Nall2000,[3 4]));
Nall2100=squeeze(mean(Nall2100,[3 4]));
N=squeeze(areaweightedmean(areaweightedmean(Nall2000,taream,2),nanmean(taream,2),1));
N1=squeeze(areaweightedmean(areaweightedmean(Nall2100,taream,2),nanmean(taream,2),1));

%% plot 12 cases surface N and % change

subplot(2,2,2)
plot(1:12,N(iorder),'ko'); hold on; plot(1:12,N1(iorder),'rd')
set(gca,'XTick',1:12); set(gca,'XTickLabels',{})
ylabel('N'); xlabel('(c)')
xlim([0.5 12.5]); 
set(gca,'fontsize',12); legend('2000','2100')
subplot(2,2,4)
plot(1:12,(N1(iorder)-N(iorder)),'k*')
hold on; plot(2,(N1(11)-N(11)),'o','Color',[0.85 0.33 0.1],'MarkerSize',9)
plot(9,(N1(4)-N(4)),'o','Color',[0 0.45 0.74],'MarkerSize',9)
ylabel('\Delta N')
set(gca,'XTick',1:12); set(gca,'XTickLabels',num2str(round(orderp2,1)))
xlabel({'case timescale','(d)'})
xlim([0.5 12.5]);
set(gca,'fontsize',12)

