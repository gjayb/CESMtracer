%making 1 figure to show climate changes in strat light and circulation

%load small
dHMXLmax=ncread('HMXLmaxdiff2111.nc','HMXL'); %using XMXL instead!
dI=ncread('SHFQSWmaxdiff2111.nc','SHF_QSW');
dSSH=0.01*ncread('sshAnnualDiff.nc','SSH');
dW=ncread('w100mAnnualDiff.nc','WVEL');

%load potential density
PD=1000*(ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','PD',[1 1 1 109],[Inf Inf 16 Inf])-1);
PD1=1000*(ncread('g.e21.G1850ECO.T62_g17.param2100.123.pop.h.PD.002101-003012.nc','PD',[1 1 1 109],[Inf Inf 16 Inf])-1);
z=ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','z_t')./100;
strat=mean(PD(:,:,16,end-11:end)-PD(:,:,6,end-11:end),4)./100;%year 10 mean of potential density at 155m-55m/100m
strat1=mean(PD1(:,:,16,end-11:end)-PD1(:,:,6,end-11:end),4)./100;
dstrat=strat1-strat;
clear PD PD1

%load XMXL, the monthly max MLD (dHMXLmax is change in annual-max
%monthly-mean MLD)
XMXL=0.01*(ncread('g.e21.G.T62_g17.param2000.121.pop.h.XMXL.002101-003012.nc','XMXL',[1 1 109],[Inf Inf Inf]));%0.01 to go cm-->m
XMXL1=0.01*(ncread('g.e21.G1850ECO.T62_g17.param2100.121.pop.h.XMXL.002101-003012.nc','XMXL',[1 1 109],[Inf Inf Inf]));
dXMXLmax=max(XMXL1,[],3)-max(XMXL,[],3);
clear XMXL XMXL1

load('globalLatlonbasin.mat')
y=sind(lat);

%% plot
figure; subplot(3,2,1); scatter(lon(:),y(:),9,dI(:),'filled');
hold on; scatter(lon(basin==0),y(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); set(gca,'XTickLabels',{});
c1=colorbar; caxis([-50 200]); cmocean('balance','pivot',0); c1.Label.String='\Delta W/m^2'; 
xlabel('(a)')

subplot(3,2,2); scatter(lon(:),y(:),9,dXMXLmax(:),'filled');
hold on; scatter(lon(basin==0),y(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); ylim([-1 1]); set(gca,'XTickLabels',{}); set(gca,'YTickLabels',{});
c2=colorbar; caxis([-250 250]); cmocean('balance','pivot',0); c2.Label.String='\Delta m'; 
xlabel('(d)')

subplot(3,2,3); scatter(lon(:),y(:),9,dSSH(:),'filled');
hold on; scatter(lon(basin==0),y(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); set(gca,'XTickLabels',{});
c3=colorbar;  cmocean('balance','pivot',0); c3.Label.String='\Delta m'; 
xlabel('(b)')

%add stratification here?
subplot(3,2,4); scatter(lon(:),y(:),9,dstrat(:)./100,'filled');
hold on; scatter(lon(basin==0),y(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); ylim([-1 1]); set(gca,'XTickLabels',{}); set(gca,'YTickLabels',{});
c4=colorbar;  cmocean('balance','pivot',0); c4.Label.String='\Delta d\sigma_\theta/dz'; 
xlabel('(e)')

subplot(3,2,5); scatter(lon(:),y(:),9,dW(:)*86400,'filled');
hold on; scatter(lon(basin==0),y(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); ylim([-1 1]); 
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlabel('longitude'); ylabel('latitude')
c5=colorbar; caxis([-35 35]); cmocean('balance','pivot',0); c5.Label.String='\Delta m/day'; 
xlabel('(c)')


