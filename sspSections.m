%figure for SSP, section through pacific, mean 208-230E, -80 to 15N
%% load
S=ncread('salt2000yr30.nc','SALT');
S1=ncread('salt2100yr30.nc','SALT');
T=ncread('temp2000yr30.nc','TEMP');
T1=ncread('temp2100yr30.nc','TEMP');

load('stratificationNutrientQ.mat', 'N*')
Na0=mean(Na0,4);
Na1=mean(Na1,4);
Nb0=mean(Nb0,4);
Nb1=mean(Nb1,4);

zm=ncread('salt2000yr30.nc','z_t')/100;
yi=1:242;
xi=221:240;
xim=230;
%% fig N 2000,2100, change
figure;
subplot(3,2,1)

%% fig T,S 2000,2100, change
figure;
subplot(3,2,1)
%% fig NTS changes only
dNa=squeeze(mean(Na1(xi,yi,1:35)-Na0(xi,yi,1:35),1));
dNb=squeeze(mean(Nb1(xi,yi,1:35)-Nb0(xi,yi,1:35),1));
dT=squeeze(mean(T1(xi,yi,1:35)-T(xi,yi,1:35),1));
dS=squeeze(mean(S1(xi,yi,1:35)-S(xi,yi,1:35),1));
figure;
subplot(2,2,1)
contourf(lat(xim,yi),-zm(1:35),dNa.',-5:0.5:5)
cmocean('balance','pivot',0)
set(gca,'fontsize',12); xlabel('(a)')
ylabel('depth (m)'); title('\Delta N slow')
xlim([-70 10])

subplot(2,2,2)
contourf(lat(xim,yi),-zm(1:35),dNb.',-5:0.5:5)
cmocean('balance','pivot',0); colorbar
set(gca,'fontsize',12); xlabel('(b)')
title('\Delta N fast')
xlim([-70 10])
subplot(2,2,3)
contourf(lat(xim,yi),-zm(1:35),dT.',-5:0.5:5)
cmocean('-tarn','pivot',0)
set(gca,'fontsize',12)
ylabel('depth (m)'); title('\Delta \theta')
xlabel({'latitude','(c)'})
xlim([-70 10])
subplot(2,2,4)
contourf(lat(xim,yi),-zm(1:35),dS.',-0.6:0.125:0.25)
cmocean('delta','pivot',0); colorbar
set(gca,'fontsize',12); xlabel('(d)')
title('\Delta Salinity')
xlim([-70 10])
