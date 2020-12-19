%Figures 3 and 4-- global N and Production spatial patterns

%load production
load('yr10JN100m.mat', 'JNall2000','JNall2100','alpha','kn','mu')
prod=squeeze(mean(sum(-JNall2000(:,:,:,:,[4 11])*(10*86400*365*1e-3*14*117/16),3),4)); %annual gC/m^2
prod1=squeeze(mean(sum(-JNall2100(:,:,:,:,[4 11])*(10*86400*365*1e-3*14*117/16),3),4));
clear JNall2*

%load regions and geometry
load('globalLatlonbasin.mat')
load('climateAndRegions.mat','osmosis','subtropSPac')
y1=sind(lat);
%% compute pattern correlation
pc1=100*(prod1(:,:,1)-prod(:,:,1))./prod(:,:,1);
pc2=100*(prod1(:,:,2)-prod(:,:,2))./prod(:,:,2);
pc1=pc1(~isnan(pc1)); pc2=pc2(~isnan(pc2));
[c1,c2]=corrcoef(pc1,pc2)

%% plot production

figure; subplot(2,2,1)
scatter(lon(:),y1(:),16,reshape(prod(:,:,1),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); cmocean('algae'); caxis([0 120]); 
set(gca,'XTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
contour(lon,y1,subtropSPac,[1 1],'m')
plot([-5 365],min(y1(basin==10)).*[1 1],'k')%arctic
contour(lon,y1,osmosis,[1 1],'c')
set(gca,'fontsize',12); xlabel('(a)'); title('slow case')

subplot(2,2,2)
scatter(lon(:),y1(:),16,reshape(prod(:,:,2),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{}); set(gca,'XTickLabels',{});
xlim([0 360]); cmocean('algae'); caxis([0 120]); 
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
contour(lon,y1,subtropSPac,[1 1],'m')
plot([-5 365],min(y1(basin==10)).*[1 1],'k')%arctic
contour(lon,y1,osmosis,[1 1],'c')
set(gca,'fontsize',12); xlabel('(b)'); title('fast case'); 
c2=colorbar; c2.Label.String='gC/m^2'; 

subplot(2,2,3)
scatter(lon(:),y1(:),16,reshape(100*(prod1(:,:,1)-prod(:,:,1))./prod(:,:,1),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); caxis([-100 200]); cmocean('delta','pivot',0)
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
contour(lon,y1,subtropSPac,[1 1],'m')
plot([-5 365],min(y1(basin==10)).*[1 1],'k')%arctic
contour(lon,y1,osmosis,[1 1],'c')
set(gca,'fontsize',12)
xlabel({'longitude','(c)'}); 

subplot(2,2,4)
scatter(lon(:),y1(:),16,reshape(100*(prod1(:,:,2)-prod(:,:,2))./prod(:,:,2),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{})
xlim([0 360]); caxis([-100 200]); cmocean('delta','pivot',0)
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
contour(lon,y1,subtropSPac,[1 1],'m')
plot([-5 365],min(y1(basin==10)).*[1 1],'k')%arctic
contour(lon,y1,osmosis,[1 1],'c')
set(gca,'fontsize',12); xlabel('(d)'); 
c4=colorbar; c4.Label.String='% change'; 
%% load nutrient
load('yr10N100m.mat', 'Nall2000','Nall2100')
N=mean(Nall2000(:,:,:,:,[4 11]),[3 4]);
N1=mean(Nall2100(:,:,:,:,[4 11]),[3 4]);
clear Nall2*
%% plot nutrients

figure; subplot(2,2,1)
scatter(lon(:),y1(:),16,reshape(N(:,:,1),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); cmocean('amp'); caxis([0 10]); 
set(gca,'XTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
contour(lon,y1,subtropSPac,[1 1],'m')
plot([-5 365],min(y1(basin==10)).*[1 1],'k')%arctic
contour(lon,y1,osmosis,[1 1],'c')
set(gca,'fontsize',12); xlabel('(a)'); title('slow case')

subplot(2,2,2)
scatter(lon(:),y1(:),16,reshape(N(:,:,2),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{}); set(gca,'XTickLabels',{});
xlim([0 360]); cmocean('amp'); caxis([0 10]); 
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
contour(lon,y1,subtropSPac,[1 1],'m')
plot([-5 365],min(y1(basin==10)).*[1 1],'k')%arctic
contour(lon,y1,osmosis,[1 1],'c')
set(gca,'fontsize',12); xlabel('(b)'); title('fast case'); 
c2=colorbar; c2.Label.String='gN'; 

subplot(2,2,3)
scatter(lon(:),y1(:),16,reshape(100*(N1(:,:,1)-N(:,:,1))./N(:,:,1),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); caxis([-100 200]); cmocean('balance','pivot',0)
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
contour(lon,y1,subtropSPac,[1 1],'m')
plot([-5 365],min(y1(basin==10)).*[1 1],'k')%arctic
contour(lon,y1,osmosis,[1 1],'c')
set(gca,'fontsize',12)
xlabel({'longitude','(c)'}); 

subplot(2,2,4)
scatter(lon(:),y1(:),16,reshape(100*(N1(:,:,2)-N(:,:,2))./N(:,:,2),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{})
xlim([0 360]); caxis([-100 200]); cmocean('balance','pivot',0)
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
contour(lon,y1,subtropSPac,[1 1],'m')
plot([-5 365],min(y1(basin==10)).*[1 1],'k')%arctic
contour(lon,y1,osmosis,[1 1],'c')
set(gca,'fontsize',12); xlabel('(d)'); 
c4=colorbar; c4.Label.String='% change'; 