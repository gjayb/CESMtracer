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

%% dprod/dN, dprod/dI
%d(muQL)dN=mu*meanL*k./(meanN+k)^2=
%d(muQL)dI=mu*meanQ*alpha*(exp(-alpha*meanI))=mu*meanQ*alpha*(1-meanL)
%for theory, assume meanQ is 0.5 
%and assume meanL is 100m-avg of L for I=100e^-z/H
H=-100./log(0.01);%so I=100*exp(-z/H) is 100 @surface, 1 @100m
z=0:100; I=100*exp(-z./H); for k=1:12; meanLtheory(k)=mean(1-exp(-alpha(k).*I)); end
dproddNtheory=mu.*meanLtheory.*(0.25./kn);
dproddItheory=0.5*mu.*alpha.*(1-meanLtheory);
prodtheory=mu.*meanLtheory.*0.5;

%try N=sqrt(kn./mu)
dproddNtheory2=mu.*meanLtheory.*kn./(kn+sqrt(kn./mu)).^2;
dproddItheory2=mu.*(sqrt(kn./mu)./(kn+sqrt(kn./mu))).*alpha.*(1-meanLtheory);
prodtheory2=mu.*meanLtheory.*(sqrt(kn./mu)./(kn+sqrt(kn./mu)));

globalmeanL=squeeze(areaweightedmean(areaweightedmean(meanL02,taream,2),mean(taream,2),1));
globalmeanQ=squeeze(areaweightedmean(areaweightedmean(meanQ02,taream,2),mean(taream,2),1));
globalmeanN=squeeze(areaweightedmean(areaweightedmean(squeeze(nanmean(Nall2000,[3 4])),taream,2),mean(taream,2),1));

globalmeanL1=squeeze(areaweightedmean(areaweightedmean(meanL12,taream,2),mean(taream,2),1));
globalmeanQ1=squeeze(areaweightedmean(areaweightedmean(meanQ12,taream,2),mean(taream,2),1));
globalmeanN1=squeeze(areaweightedmean(areaweightedmean(squeeze(nanmean(Nall2100,[3 4])),taream,2),mean(taream,2),1));


dproddN=mu(:).*globalmeanL.*kn(:)./((globalmeanN+kn(:)).^2);
dproddI=mu(:).*globalmeanQ.*alpha(:).*(1-globalmeanL);
%% plot
figure; subplot(2,2,1); semilogx(dproddNtheory,1:12,'o'); hold on; semilogx(dproddItheory,1:12,'o')
ylabel('case'); legend('dR/dN Q=0.5','dR/dI N=k')
subplot(2,2,2); semilogx(dproddNtheory2,1:12,'o'); hold on; semilogx(dproddItheory2,1:12,'o')
legend('dR/dN N=(k/mu)^{0.5}','dR/dI')
subplot(2,2,3); semilogx(dproddNtheory,100*(prod1-prod)./prod,'o'); hold on; semilogx(dproddItheory,100*(prod1-prod)./prod,'o');
ylabel('% change global R')
subplot(2,2,4); semilogx(dproddNtheory2,100*(prod1-prod)./prod,'o'); hold on; semilogx(dproddItheory2,100*(prod1-prod)./prod,'o');

figure; subplot(2,2,1); semilogx(dproddN,1:12,'o'); hold on; semilogx(dproddI,1:12,'o')
subplot(2,2,2); loglog(dproddN,dproddNtheory,'o'); hold on; loglog(dproddI,dproddItheory,'o')
loglog([3e-4 1],[3e-4 1],'k')
subplot(2,2,4)
loglog(dproddN,dproddNtheory2,'o'); hold on; loglog(dproddI,dproddItheory2,'o')
loglog([3e-4 1],[3e-4 1],'k')
subplot(2,2,3)
semilogx(dproddN,100*(prod1-prod)./prod,'o'); hold on; semilogx(dproddI,100*(prod1-prod)./prod,'o');
%% compute I from L?
L02=L0; L02(L02>1)=1; L02(L02<0)=0;
alpha5(1,1,1,1,1:12)=alpha;
I0=-log(1-L02)./repmat(alpha5,[320 384 10 12 1]);
I0(isinf(I0))=NaN;

L12=L0; L12(L12>1)=1; L12(L12<0)=0;
I1=-log(1-L12)./repmat(alpha5,[320 384 10 12 1]);

meanI0=squeeze(nanmean(I0,3));
meanI1=squeeze(nanmean(I1,3));
meanN0=squeeze(nanmean(Nall2000,3));
meanN1=squeeze(nanmean(Nall2100,3));
%% frequency of N and I, production rate %NOT VERY HELPFUL

edgeI=[0:0.01:0.1 0.2:0.1:1 2:10 20:10:100]; % 0.2:5 5.5:0.5:20 21:65];
edgeN=[0:0.01:0.1 0.2:0.1:1 2:20];%[0:0.1:4 5:20];
binI=0.5*(edgeI(2:end)+edgeI(1:end-1));
binN=0.5*(edgeN(2:end)+edgeN(1:end-1));
[gN,gI]=meshgrid(binN,binI);
%
%freqNI=NaN(length(binI),length(binN),12);
%freqNI1=NaN(length(binI),length(binN),12);

freqNImean=NaN(length(binI),length(binN),12);
freqNImean1=NaN(length(binI),length(binN),12);

for k=[4 11]
    %holdI0=reshape(I0(:,:,:,:,k),1,[]);
    %holdN0=reshape(Nall2000(:,:,:,:,k),1,[]);
    %holdI1=reshape(I1(:,:,:,:,k),1,[]);
    %holdN1=reshape(Nall2100(:,:,:,:,k),1,[]);
        holdmeanI0=reshape(meanI0(:,:,:,k),1,[]);
    holdmeanN0=reshape(meanN0(:,:,:,k),1,[]);
    holdmeanI1=reshape(meanI1(:,:,:,k),1,[]);
    holdmeanN1=reshape(meanN1(:,:,:,k),1,[]);
    for i=1:length(binI)
        i
        for j=1:length(binN)
            %freqNI(i,j,k)=sum(holdI0>=edgeI(i)&holdI0<edgeI(i+1)&holdN0>=edgeN(j)&holdN0<edgeN(j+1));
            %freqNI1(i,j,k)=sum(holdI1>=edgeI(i)&holdI1<edgeI(i+1)&holdN1>=edgeN(j)&holdN1<edgeN(j+1));
            freqNImean(i,j,k)=sum(holdmeanI0>=edgeI(i)&holdmeanI0<edgeI(i+1)&holdmeanN0>=edgeN(j)&holdmeanN0<edgeN(j+1));
            freqNImean1(i,j,k)=sum(holdmeanI1>=edgeI(i)&holdmeanI1<edgeI(i+1)&holdmeanN1>=edgeN(j)&holdmeanN1<edgeN(j+1));
        end
    end
end
%
clear totNI* prodNI
totNI=nansum(nansum(freqNI));
totNI1=nansum(nansum(freqNI1));


for k=1:12
prodNI(:,:,k)=mu(k)*(gN./(gN+kn(k))).*(1-exp(-alpha(k)*gI));
justQL(:,:,k)=(gN./(gN+kn(k))).*(1-exp(-alpha(k)*gI));
end

%% plot N,I, prod, freq %NOT VERY HELPFUL

figure; subplot(1,2,1)
pcolor(binN,binI,justQL(:,:,4)); shading 'flat'; 
cmocean('haline'); c1=colorbar; c1.Label.String='gN/day';
hold on
[cca,ccb]=contour(binN,binI,freqNI(:,:,4)/totNI(4),'k'); clabel(cca,ccb);
[ccc,ccd]=contour(binN,binI,freqNI1(:,:,4)/totNI1(4),'w'); clabel(ccc,ccd);
title('Slow case'); xlabel('N'); ylabel('I')
set(gca,'YScale','log'); set(gca,'XScale','log')

subplot(1,2,2)
pcolor(binN,binI,justQL(:,:,11)); shading 'flat'; 
cmocean('haline'); c1=colorbar; c1.Label.String='gN/day';
hold on
[cca,ccb]=contour(binN,binI,freqNI(:,:,11)/totNI(11),'k'); clabel(cca,ccb);
[ccc,ccd]=contour(binN,binI,freqNI1(:,:,11)/totNI1(11),'w'); clabel(ccc,ccd);
title('Fast case'); xlabel('N'); ylabel('I')
set(gca,'YScale','log'); set(gca,'XScale','log')
%% compute Q, L frequency
edgeQ=0:0.05:1;
edgeL=0:0.05:1;
binQ=0.5*(edgeQ(1:end-1)+edgeQ(2:end));
binL=0.5*(edgeL(1:end-1)+edgeL(2:end));
[Qg,Lg]=meshgrid(binQ,binL);

freqQL=NaN(length(binQ),length(binL),12);
freqQL1=freqQL; freqQLmean=freqQL; freqQLmean1=freqQL;

for k=1:12
    k
    holdQ0=reshape(Q0(:,:,:,:,k),1,[]);
    holdL0=reshape(L0(:,:,:,:,k),1,[]);
    holdQ1=reshape(Q1(:,:,:,:,k),1,[]);
    holdL1=reshape(L1(:,:,:,:,k),1,[]);
    holdmeanQ0=reshape(meanQ0(:,:,:,k),1,[]);
    holdmeanL0=reshape(meanL0(:,:,:,k),1,[]);
    holdmeanQ1=reshape(meanQ1(:,:,:,k),1,[]);
    holdmeanL1=reshape(meanL1(:,:,:,k),1,[]);
    for i=1:length(binQ)
        %i
        for j=1:length(binL)
            freqQL(i,j,k)=sum(holdQ0>=edgeQ(i)&holdQ0<edgeQ(i+1)&holdL0>=edgeL(j)&holdL0<edgeL(j+1));
            freqQL1(i,j,k)=sum(holdQ1>=edgeQ(i)&holdQ1<edgeQ(i+1)&holdL1>=edgeL(j)&holdL1<edgeL(j+1));
            freqQLmean(i,j,k)=sum(holdmeanQ0>=edgeQ(i)&holdmeanQ0<edgeQ(i+1)&holdmeanL0>=edgeL(j)&holdmeanL0<edgeL(j+1));
            freqQLmean1(i,j,k)=sum(holdmeanQ1>=edgeQ(i)&holdmeanQ1<edgeQ(i+1)&holdmeanL1>=edgeL(j)&holdmeanL1<edgeL(j+1));
        end
    end
end
tot1=nansum(nansum(freqQL(:,:,1)));
tot2=nansum(nansum(freqQLmean(:,:,1)));
%% plot Q,L frequency %depth-avg ones show more info
figure; subplot(2,2,1)
pcolor(Qg,Lg,freqQL(:,:,4)/tot1); shading 'flat'; hold on
[cc1,cc2]=contour(Qg,Lg,freqQL1(:,:,4)/tot1,'w'); clabel(cc1,cc2);
title('slow'); ylabel({'each depth','L'}); 

subplot(2,2,2)
pcolor(Qg,Lg,freqQL(:,:,11)/tot1); shading 'flat'; hold on
[cc1,cc2]=contour(Qg,Lg,freqQL1(:,:,11)/tot1,'w'); clabel(cc1,cc2);
title('fast')

subplot(2,2,3)
pcolor(Qg,Lg,freqQLmean(:,:,4)/tot2); shading 'flat'; hold on
[cc1,cc2]=contour(Qg,Lg,freqQLmean1(:,:,4)/tot2,'w'); clabel(cc1,cc2);
ylabel({'100m depth-avg','L'}); xlabel('Q')

subplot(2,2,4)
pcolor(Qg,Lg,freqQLmean(:,:,11)/tot2); shading 'flat'; hold on
[cc1,cc2]=contour(Qg,Lg,freqQLmean1(:,:,11)/tot2,'w'); clabel(cc1,cc2);

%% compute deltaQL=deltaQ*L+deltaL*Q+deltaQ*deltaL+residual, Q1/Q0, L1/L0

dQL=squeeze(mean(meanQ1.*meanL1-meanQ0.*meanL0,3));
LdQ=squeeze(mean(meanL0.*(meanQ1-meanQ0),3));
QdL=squeeze(mean(meanQ0.*(meanL1-meanL0),3));
dQdL=squeeze(mean((meanL1-meanL0).*(meanQ1-meanQ0),3));
reynoldsresid=dQL-LdQ-QdL-dQdL;%zero to numerical precision (1e-16)
Qratio=meanQ12./meanQ02;
Lratio=meanL12./meanL02;
QL0=squeeze(mean(meanQ0.*meanL0,3));

%correlations of components with dQL
for k=1:12
    v1=dQL(:,:,k); x=v1(~isnan(v1));
    y=LdQ(:,:,k); y=y(~isnan(y));
   hold1=corrcoef(x,y);
   cordQL(1,k)=hold1(2); %1 is LdQ, 2 is QdL, 3 is dQdL
    y=QdL(:,:,k); y=y(~isnan(y));
   hold1=corrcoef(x,y);
   cordQL(2,k)=hold1(2);
   y=dQdL(:,:,k); y=y(~isnan(y));
   hold1=corrcoef(x,y);
   cordQL(3,k)=hold1(2);
end

%% correlations of components across parameter space
for k=1:12
    corAll(k,k,1:4)=1;
    for i=k+1:12
    x=dQL(:,:,k); x=x(~isnan(x));
    y=dQL(:,:,i); y=y(~isnan(y));
   hold1=corrcoef(x,y);
   corAll(k,i,1)=hold1(2); %1 is dQL 2 is LdQ, 3 is QdL, 4 is dQdL
   corAll(i,k,1)=hold1(2);
    x=LdQ(:,:,k); x=x(~isnan(x));
    y=LdQ(:,:,i); y=y(~isnan(y));
   hold1=corrcoef(x,y);
   corAll(k,i,2)=hold1(2);
   corAll(i,k,2)=hold1(2);
   x=QdL(:,:,k); x=x(~isnan(x));
    y=QdL(:,:,i); y=y(~isnan(y));
   hold1=corrcoef(x,y);
   corAll(k,i,3)=hold1(2);
   corAll(i,k,3)=hold1(2);
    x=dQdL(:,:,k); x=x(~isnan(x));
    y=dQdL(:,:,i); y=y(~isnan(y));
    hold1=corrcoef(x,y);
   corAll(k,i,4)=hold1(2);
   corAll(i,k,4)=hold1(2);
    end
end
%
for k=1:12
    corQL(k,k,1:4)=1;
    for i=k+1:12
    x=meanQ02(:,:,k); x=x(~isnan(x));
    y=meanQ02(:,:,i); y=y(~isnan(y));
   hold1=corrcoef(x,y);
   corQL(k,i,1)=hold1(2); %1 is Q, 2 is L, 3 is dQ, 4 is dL
   corQL(i,k,1)=hold1(2);
    x=meanL02(:,:,k); x=x(~isnan(x));
    y=meanL02(:,:,i); y=y(~isnan(y));
   hold1=corrcoef(x,y);
   corQL(k,i,2)=hold1(2);
   corQL(i,k,2)=hold1(2);
   x=meanQ12(:,:,k)-meanQ02(:,:,k); x=x(~isnan(x));
    y=meanQ12(:,:,i)-meanQ02(:,:,i); y=y(~isnan(y));
   hold1=corrcoef(x,y);
   corQL(k,i,3)=hold1(2);
   corQL(i,k,3)=hold1(2);
    x=meanL12(:,:,k)-meanL02(:,:,k); x=x(~isnan(x));
    y=meanL12(:,:,i)-meanL02(:,:,i); y=y(~isnan(y));
    hold1=corrcoef(x,y);
   corQL(k,i,4)=hold1(2);
   corQL(i,k,4)=hold1(2);
    end
end

%% plot QL, delta(QL), deltaQ*L, deltaL*Q, deltaQ*deltaL 
figure; subplot(3,2,1) %slow
scatter(lon(:),y1(:),16,reshape(QL0(:,:,4),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{'-60','-30','0','30','60'})
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); %caxis([0 1]); 
cmocean('algae'); ylabel('QL_{2000}'); title('slow case'); colorbar

subplot(3,2,2)
scatter(lon(:),y1(:),16,reshape(dQL(:,:,4),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); caxis([-0.1 0.05]); cmocean('delta','pivot',0); colorbar
cmocean('delta','pivot',0); ylabel('QL_{2100}-QL_{2000}'); 

subplot(3,2,3)
scatter(lon(:),y1(:),16,reshape(LdQ(:,:,4),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{'-60','-30','0','30','60'})
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); caxis([-0.1 0.05]); cmocean('delta','pivot',0); colorbar
cmocean('delta','pivot',0); ylabel('L_{2000}*(Q_{2100}-Q_{2000})'); 

subplot(3,2,4)
scatter(lon(:),y1(:),16,reshape(QdL(:,:,4),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); caxis([-0.1 0.05]); cmocean('delta','pivot',0); colorbar
cmocean('delta','pivot',0); ylabel('Q_{2000}*(L_{2100}-L_{2000})'); 

subplot(3,2,5)
scatter(lon(:),y1(:),16,reshape(dQdL(:,:,4),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{'-60','-30','0','30','60'})
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); caxis([-0.1 0.05]); cmocean('delta','pivot',0); colorbar
cmocean('delta','pivot',0); ylabel('\DeltaQ*\DeltaL'); 

subplot(3,2,6)
scatter(lon(:),y1(:),16,reshape(reynoldsresid(:,:,4),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); colorbar %caxis([-1 1]); 
cmocean('delta','pivot',0); ylabel('residual'); 
%
figure; subplot(3,2,1) %fast
scatter(lon(:),y1(:),16,reshape(QL0(:,:,11),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{'-60','-30','0','30','60'})
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); %caxis([0 1]); 
cmocean('algae'); ylabel('QL_{2000}'); title('fast case'); colorbar

subplot(3,2,2)
scatter(lon(:),y1(:),16,reshape(dQL(:,:,11),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); %caxis([-0.1 0.05]); 
cmocean('delta','pivot',0); colorbar
 ylabel('QL_{2100}-QL_{2000}'); 

subplot(3,2,3)
scatter(lon(:),y1(:),16,reshape(LdQ(:,:,11),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{'-60','-30','0','30','60'})
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); %caxis([-0.1 0.05]); 
cmocean('delta','pivot',0); colorbar
 ylabel('L_{2000}*(Q_{2100}-Q_{2000})'); 

subplot(3,2,4)
scatter(lon(:),y1(:),16,reshape(QdL(:,:,11),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); %caxis([-0.1 0.05]); 
cmocean('delta','pivot',0); colorbar
 ylabel('Q_{2000}*(L_{2100}-L_{2000})'); 

subplot(3,2,5)
scatter(lon(:),y1(:),16,reshape(dQdL(:,:,11),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{'-60','-30','0','30','60'})
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); %caxis([-0.1 0.05]); 
cmocean('delta','pivot',0); colorbar
 ylabel('\DeltaQ*\DeltaL'); 

subplot(3,2,6)
scatter(lon(:),y1(:),16,reshape(reynoldsresid(:,:,11),1,[]),'filled')
set(gca,'YTick',[-0.866 -0.5 0 0.5 0.866])
set(gca,'YTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
xlim([0 360]); colorbar %caxis([-1 1]); 
cmocean('delta','pivot',0); ylabel('residual'); 

%% plot mean Q, mean L, Q1/Q0, L1/L0 
figure; subplot(2,2,1)%mean Q,L, Q1/Q0,L1/L0, slow case
scatter(lon(:),y1(:),16,reshape(meanQ02(:,:,4),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); caxis([0 1]); 
set(gca,'XTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); xlabel('(a)'); title('slow case Q')

subplot(2,2,2)
scatter(lon(:),y1(:),16,reshape(meanL02(:,:,4),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{}); set(gca,'XTickLabels',{});
xlim([0 360]);  caxis([0 1]);
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); xlabel('(b)'); title('slow case L'); 
c2=colorbar; %c2.Label.String='months Q<L'; 

subplot(2,2,3)
scatter(lon(:),y1(:),16,reshape((Qratio(:,:,4)),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); %caxis([-12.5 12.5]); 
%cmocean('-delta',25,'pivot',0)
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); ylabel('latitude');
xlabel({'longitude','(c)'});  title('Q_{2100}/Q_{2000}')

subplot(2,2,4)
scatter(lon(:),y1(:),16,reshape((Lratio(:,:,4)),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{})
xlim([0 360]); %caxis([-12.5 12.5]); 
%cmocean('-delta',25,'pivot',0)
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); xlabel('(d)'); title('L_{2100}/L_{2000}')
c4=colorbar; %c4.Label.String='\Delta months'; 

%
figure; subplot(2,2,1)%mean Q,L, Q1/Q0,L1/L0, slow case
scatter(lon(:),y1(:),16,reshape(meanQ02(:,:,11),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); caxis([0 1]); 
set(gca,'XTickLabels',{});
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); xlabel('(a)'); title('fast case Q')

subplot(2,2,2)
scatter(lon(:),y1(:),16,reshape(meanL02(:,:,11),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{}); set(gca,'XTickLabels',{});
xlim([0 360]);  caxis([0 1]);
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); xlabel('(b)'); title('fast case L'); 
c2=colorbar; %c2.Label.String='months Q<L'; 

subplot(2,2,3)
scatter(lon(:),y1(:),16,reshape((Qratio(:,:,11)),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{'-90','-60','-30','0','30','60','90'})
xlim([0 360]); %caxis([-12.5 12.5]); 
%cmocean('-delta',25,'pivot',0)
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); ylabel('latitude');
xlabel({'longitude','(c)'});  title('Q_{2100}/Q_{2000}')

subplot(2,2,4)
scatter(lon(:),y1(:),16,reshape((Lratio(:,:,11)),1,[]),'filled')
set(gca,'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1])
set(gca,'YTickLabels',{})
xlim([0 360]); %caxis([-12.5 12.5]); 
%cmocean('-delta',25,'pivot',0)
hold on; scatter(lon(basin==0),y1(basin==0),4,[0.7 0.7 0.7],'filled')
set(gca,'Color',[0.7 0.7 0.7]);
set(gca,'fontsize',12); xlabel('(d)'); title('L_{2100}/L_{2000}')
c4=colorbar; %c4.Label.String='\Delta months'; 
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
dPdL(i)=p(1);
ndPdL(i)=length(x);
end
%% compute &plot semitheoretical dP/dQ and dP/dL
%dP/dQ=mu*L %dP/dL=mu*Q 
%note mu in days already, doing rest of standard unit conversion 
   dPdQtheory=(10*365*1e-3*14*117/16)*mu(:).*squeeze(areaweightedmean(areaweightedmean(meanL02,taream,2),mean(taream,2),1));
   dPdLtheory=(10*365*1e-3*14*117/16)*mu(:).*squeeze(areaweightedmean(areaweightedmean(meanQ02,taream,2),mean(taream,2),1));    

figure; 
yyaxis right
plot(dPdLtheory); hold on; scatter(1:12,dPdLtheory,16,globalmeanQltL0,'filled')
caxis([-0.5 12.5]); c3=colorbar; c3.Label.String='months Q<L';
ylabel('dProduction/dL=mu Q')
yyaxis left
plot(dPdQtheory); hold on; scatter(1:12,dPdQtheory,16,globalmeanQltL0,'filled')
colormap(cmap1(end:-1:1,:)); caxis([-0.5 12.5]); 
ylabel('dProduction/dQ=mu L'); xlabel('parameter cases')
xlim([0.5 12.5])
set(gca,'fontsize',12)

%% fig 6, 
% a)slow case delta P vs delta L color delta Q
% b)fast case delta P vs delta Q color delta L
% c)x 12 cases y1 dP/dQ y2 dP/dL color months Q<L
cmap1=parula(13);
f1=figure(2);

subplot(3,1,3)
yyaxis right
plot(dPdLb); hold on; scatter(1:12,dPdLb,16,globalmeanQltL0,'filled')
plot(dPdLtheory,'--o')
caxis([-0.5 12.5]); c3=colorbar; c3.Label.String='months Q<L';
ylabel('dProduction/dL')
yyaxis left
plot(dPdQb); hold on; scatter(1:12,dPdQb,16,globalmeanQltL0,'filled')
plot(dPdQtheory,'--o')
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

