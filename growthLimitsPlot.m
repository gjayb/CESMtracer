% paper figs 5 and 6 on growth limits
% fig 5, spatial patterns of limiting factor
% fig 6, global delta P vs delta L or delta Q examples and 12-case dP/dQ,dPdL

%load
load('yr10JN100m.mat')
load('yr10N100m.mat','Nall*')
load('globalLatlonbasin.mat')
load('choosingParamsNnew.mat', 'tarea')
taream=tarea*1e-4;
kn5(1,1,1,1,:)=kn;
mu5(1,1,1,1,:)=mu./86400;

Q0=Nall2000./(Nall2000+repmat(kn5,[320 384 10 12 1]));
Q1=Nall2100./(Nall2100+repmat(kn5,[320 384 10 12 1]));
L0=-JNall2000./(repmat(mu5,[320 384 10 12 1]).*Q0);
L1=-JNall2100./(repmat(mu5,[320 384 10 12 1]).*Q1);

meanQ0=squeeze(nanmean(Q0,3));
meanQ1=squeeze(nanmean(Q1,3));
meanL0=squeeze(nanmean(L0,3));
meanL1=squeeze(nanmean(L1,3));
QltL0=squeeze(sum(meanQ0<meanL0,3));
QltL1=squeeze(sum(meanQ1<meanL1,3));
QltL1(repmat(basin,[1 1 12])==0)=NaN;
QltL0(repmat(basin,[1 1 12])==0)=NaN;

meanL02=squeeze(nanmean(meanL0,3));
meanL12=squeeze(nanmean(meanL1,3));
meanQ02=squeeze(nanmean(meanQ0,3));
meanQ12=squeeze(nanmean(meanQ1,3));

globalmeanQltL0=squeeze(areaweightedmean(areaweightedmean(QltL0,taream,2),mean(taream,2),1));
globalmeanQltL1=squeeze(areaweightedmean(areaweightedmean(QltL1,taream,2),mean(taream,2),1));

meanJN0=squeeze(mean(sum((-10*86400*365*1e-3*14*117/16)*JNall2000,3),4)); %10*86400*365*1e-3*14*117/16
meanJN1=squeeze(mean(sum((-10*86400*365*1e-3*14*117/16)*JNall2100,3),4));

y1=sind(lat);
%% Plot spatial patterns of months limited
cmap1=parula(13);
close all
figure(1)
subplot(2,2,1); 
scatter(lon(:),y1(:),16,reshape(QltL0(:,:,4),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); colormap(cmap1(end:-1:1,:)); caxis([-0.5 12.5]); 
set(gca,'XTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); xlabel('(a)'); title('slow case')

subplot(2,2,2)
scatter(lon(:),y1(:),16,reshape(QltL0(:,:,11),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{}); set(gca,'XTickLabels',{});
xlim([0 360]); colormap(parula(13)); caxis([-0.5 12.5]);
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); xlabel('(b)'); title('fast case'); 
c2=colorbar; c2.Label.String='months Q<L'; 

subplot(2,2,3)
scatter(lon(:),y1(:),16,reshape(QltL1(:,:,4)-QltL0(:,:,4),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); caxis([-12.5 12.5]); cmocean('-delta',25,'pivot',0)
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12)
xlabel({'longitude','(c)'}); 

subplot(2,2,4)
scatter(lon(:),y1(:),16,reshape(QltL1(:,:,11)-QltL0(:,:,11),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{})
xlim([0 360]); caxis([-12.5 12.5]); cmocean('-delta',25,'pivot',0)
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); xlabel('(d)'); 
c4=colorbar; c4.Label.String='\Delta months'; 

%% compute dP/dL, dQ/dL all locations
for i=1:12
    
holdvar1=meanJN1(:,:,(i))-meanJN0(:,:,(i));
holdvar2=meanQ12(:,:,(i))-meanQ02(:,:,(i));
holdvar3=meanL12(:,:,(i))-meanL02(:,:,(i));

holdlogic2=abs(holdvar3)<0.005;%change in L is small
    x=holdvar2(holdlogic2&~isnan(holdvar1)); y=holdvar1(holdlogic2&~isnan(holdvar1));
p=polyfit(x,y,1);
f=polyval(p,-0.5:0.05:0.5);
figure; plot(x,y,'.')
hold all; plot(-0.5:0.05:0.5,f,'k')
title(num2str(i)); xlabel('dQ')
dPdQb(i)=p(1);
ndPdQb(i)=length(x);

holdlogic2=abs(holdvar2)<0.005;%change in Q is small
x=holdvar3(holdlogic2&~isnan(holdvar1)); y=holdvar1(holdlogic2&~isnan(holdvar1));
p=polyfit(x,y,1);
f=polyval(p,-0.5:0.05:0.5);
figure; plot(x,y,'.')
hold all; plot(-0.5:0.05:0.5,f,'k')
title(num2str(i)); xlabel('dL')
dPdLb(i)=p(1);
ndPdLb(i)=length(x);
end
%% compute dP/dQ, dP/dL from Q,L-limited regions
%this is the one used in the plot!!
for i=1:12
    holdlogic=(QltL0(:,:,i)>9);%&(QltL1(:,:,i)>9);%Q<L
holdvar1=meanJN1(:,:,(i))-meanJN0(:,:,(i));
holdvar2=meanQ12(:,:,(i))-meanQ02(:,:,(i));
holdvar3=meanL12(:,:,(i))-meanL02(:,:,(i));

holdlogic2=abs(holdvar3)<0.01;%change in L is small
    x=holdvar2(holdlogic&holdlogic2&~isnan(holdvar1)); y=holdvar1(holdlogic&holdlogic2&~isnan(holdvar1));
p=polyfit(x,y,1);
%f=polyval(p,-0.5:0.05:0.5);
%figure; plot(x,y,'.')
%hold all; plot(-0.5:0.05:0.5,f,'k')
dPdQ(i)=p(1);
ndPdQ(i)=length(x);

holdlogic=QltL0(:,:,i)<3;%Q>L

holdlogic2=abs(holdvar2)<0.01;%change in Q is small
    x=holdvar3(holdlogic&holdlogic2&~isnan(holdvar1)); y=holdvar1(holdlogic&holdlogic2&~isnan(holdvar1));
p=polyfit(x,y,1);
f=polyval(p,-0.5:0.05:0.5);
%figure; plot(x,y,'.')
%hold all; plot(-0.5:0.05:0.5,f,'k')
%dPdL(i)=p(1);
ndPdL(i)=length(x);
end
%% fig 6, 
% a)slow case delta P vs delta L color delta Q
% b)fast case delta P vs delta Q color delta L
% c)x 12 cases y1 dP/dQ y2 dP/dL color months Q<L

f1=figure(2);

subplot(3,1,3)
yyaxis right
plot(dPdLb); hold on; scatter(1:12,dPdLb,16,globalmeanQltL0,'filled')
caxis([-0.5 12.5]); c3=colorbar; c3.Label.String='months Q<L';
ylabel('dProduction/dL')
yyaxis left
plot(dPdQb); hold on; scatter(1:12,dPdQb,16,globalmeanQltL0,'filled')
colormap(cmap1(end:-1:1,:)); caxis([-0.5 12.5]); 
ylabel('dProduction/dQ'); xlabel('parameter cases')
xlim([0.5 12.5])
set(gca,'fontsize',12)

subplot(3,1,1)
holdvar1=meanJN1(:,:,(4))-meanJN0(:,:,(4));
holdvar2=meanQ12(:,:,(4))-meanQ02(:,:,(4));
holdvar3=meanL12(:,:,(4))-meanL02(:,:,(4));
holdlogic=(QltL0(:,:,4)<6);
scatter(holdvar3(holdlogic),holdvar1(holdlogic),9,holdvar2(holdlogic),'filled')
holdlogic2=abs(holdvar2)<0.005;%change in Q is small
x=holdvar3(holdlogic&holdlogic2&~isnan(holdvar1)); y=holdvar1(holdlogic&holdlogic2&~isnan(holdvar1));
p=polyfit(x,y,1);
f=polyval(p,-0.07:0.001:0.05);
hold all; plot(-0.07:0.001:0.05,f,'k')
grid on
set(gca,'Color',[0.7 0.7 0.7]); set(gca,'fontsize',12)
cmocean('balance','pivot',0); c1=colorbar; c1.Label.String='\Delta Q';
ylabel('\Delta Production, gC/m^2'); xlabel('\Delta L')

subplot(3,1,2)
holdvar1=meanJN1(:,:,(11))-meanJN0(:,:,(11));
holdvar2=meanQ12(:,:,(11))-meanQ02(:,:,(11));
holdvar3=meanL12(:,:,(11))-meanL02(:,:,(11));
holdlogic=(QltL0(:,:,11)>9);
scatter(holdvar2(holdlogic),holdvar1(holdlogic),9,holdvar3(holdlogic),'filled')
holdlogic2=abs(holdvar3)<0.005;%change in L is small
x=holdvar2(holdlogic&holdlogic2&~isnan(holdvar1)); y=holdvar1(holdlogic&holdlogic2&~isnan(holdvar1));
p=polyfit(x,y,1);
f=polyval(p,-0.2:0.01:0.25);
hold all; plot(-0.2:0.01:0.25,f,'k')
grid on
set(gca,'Color',[0.7 0.7 0.7]); set(gca,'fontsize',12)
cmocean('delta','pivot',0); c2=colorbar; c2.Label.String='\Delta L';
xlabel('\Delta Q')

