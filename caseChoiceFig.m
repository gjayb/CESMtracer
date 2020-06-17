%Figure 2:
%annual zonal surface mean plots of N, P for 27 cases
%annual total new production in both climates for 12 cases

%load zonal annual means
load('choosingParamsNnew.mat', 'zonalNO3woa','zonalN','zonalNb','latMean')
load('choosingParamsPnew.mat','zonalP','zonalBEC_C')
zonalPb=reshape(zonalP,[58 60 27]);
better12=zeros([3 3 3]);
better12(1,2,:)=1; better12(2,1,:)=1; better12(3,1,:)=1; better12(3,3,:)=1;
%better12 [2:4 9 11:13 18 20:22 27]
slowcase=zeros([3 3 3]); slowcase(1,2,1)=1; %case 4
fastcase=zeros([3 3 3]); fastcase(3,3,2)=1; %case 18

%% plot zonal annual means
dark1=logical(reshape(better12,27,1));
figure; subplot(2,3,1); 
plot(zonalNO3woa(:,1),latMean,'k','Linewidth',2); hold on
plot(zonalNb(:,1,4),latMean,'b','Linewidth',2); plot(zonalNb(:,1,18),latMean,'r','Linewidth',2);
plot(squeeze(zonalNb(:,1,[1:3 5:17 19:27])),latMean,'Color',[0.6 0.6 0.6]); 
plot(zonalNO3woa(:,1),latMean,'k','Linewidth',2); hold on
plot(zonalNb(:,1,4),latMean,'b','Linewidth',2); plot(zonalNb(:,1,18),latMean,'r','Linewidth',2);
legend('WOA NO_3','slow','fast','25 other cases'); ylim([-80 90]); set(gca,'fontsize',12)
ylabel('latitude'); xlabel({'N units','(a)'})

subplot(2,3,2) %redfield in BEC is 16N:117C
plot(zonalBEC_C(:,1).*16./117,latMean,'k','Linewidth',2); hold on
plot(zonalPb(:,1,4),latMean,'b','Linewidth',2); plot(zonalPb(:,1,18),latMean,'r','Linewidth',2);
plot(squeeze(zonalPb(:,1,2)),latMean,'Color',[0.3 0.3 0.3]); 
plot(squeeze(zonalPb(:,1,1)),latMean,'Color',[0.8 0.8 0.8]); 
plot(squeeze(zonalPb(:,1,[5:8 10 14:17 19 23:26])),latMean,'Color',[0.8 0.8 0.8]); 
plot(squeeze(zonalPb(:,1,[3:4 9 11:13 18 20:22 27])),latMean,'Color',[0.3 0.3 0.3]); 
plot(zonalBEC_C(:,1).*16./117,latMean,'k','Linewidth',2); hold on
plot(zonalPb(:,1,4),latMean,'b','Linewidth',2); plot(zonalPb(:,1,18),latMean,'r','Linewidth',2);
legend('BEC','slow','fast','10 next best cases','15 other cases'); ylim([-80 90]); set(gca,'fontsize',12)
xlabel({'N units','(b)'}); set(gca,'YTickLabels',{})
%% load 12 cases production
load('yr10JN100m.mat', 'JNall2000','JNall2100','alpha','kn','mu')
taream=ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','TAREA')./1e4;%square cm to square m
%Jnutri100day=squeeze(sum(repmat(dz2,[320 384 1 120]).*repmat(taream,[1 1 10 120]).*Jnutri(:,:,1:10,:)*86400*365*1e-3*14,3));
%J_NUTRI is the rate mmol/(s m^3) of transfer N-->P
%convert to PgC/yr: J_NUTRI*(dz*taream m^3)*(86400*365 s/yr)*(1e-3mol/mmol)*(14gN/mol)*(117C/16N)
prod=squeeze(mean(nansum(-JNall2000*(10*86400*365*1e-3*14*117/16).*repmat(taream,[1 1 10 12 12]),[3 1 2]),4));
prod1=squeeze(mean(nansum(-JNall2100*(10*86400*365*1e-3*14*117/16).*repmat(taream,[1 1 10 12 12]),[3 1 2]),4));
clear JNall2*
%% plot 12 cases production and % change
subplot(2,3,3)
plot(1:12,prod,'o'); hold on; plot(1:12,prod1,'d')
set(gca,'XTick',1:12); set(gca,'XTickLabels',{})
ylabel('new production, gC')
xlim([0.5 12.5]); ylim([4.5e15 8e15])
set(gca,'fontsize',12)
subplot(2,3,6)
plot(1:12,100*(prod1-prod)./prod,'*')
ylabel('% change')
set(gca,'XTick',1:12); set(gca,'XTickLabels',{'1','2','3','slow','5','6','7','8','9','10','fast','12'})
xlabel({'12 cases','(c)'})
xlim([0.5 12.5]); ylim([-20 0])
set(gca,'fontsize',12)
