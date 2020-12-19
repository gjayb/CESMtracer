%%
load('yr10JN100m.mat', 'JNall2000','JNall2100','alpha','kn','mu')
taream=ncread('g.e21.G.T62_g17.param2000.123.pop.h.PD.002101-003012.nc','TAREA')./1e4;%square cm to square m
%Jnutri100day=squeeze(sum(repmat(dz2,[320 384 1 120]).*repmat(taream,[1 1 10 120]).*Jnutri(:,:,1:10,:)*86400*365*1e-3*14,3));
%J_NUTRI is the rate mmol/(s m^3) of transfer N-->P
%convert to PgC/yr: J_NUTRI*(dz*taream m^3)*(86400*365 s/yr)*(1e-3mol/mmol)*(14gN/mol)*(117C/16N)
prod=squeeze(nansum(-JNall2000*(10*86400*365*1e-3*14*117/16).*repmat(taream,[1 1 10 12 12]),3));%xyz time case
prod1=squeeze(nansum(-JNall2100*(10*86400*365*1e-3*14*117/16).*repmat(taream,[1 1 10 12 12]),3));
clear JNall2*

orderparam12=kn(:)./mu(:)+1./(alpha(:).*mu(:));
[orderp2,iorder]=sort(orderparam12);
orderp3=orderp2; orderp3(9)=orderp3(9)+1;  orderp3(11)=orderp3(11)+1;

%%
prodm=squeeze(mean(prod,3));
prod1m=squeeze(mean(prod1,3));
load('climateAndRegions.mat', 'subtropSPac','osmosis')%,'arctic')
load('globalLatlonbasin.mat')
arctic=(lat>66.5)&(basin~=0);
for i=1:12; holdvar=prodm(:,:,i); prodtot(i,1)=nansum(holdvar(:)); end
for i=1:12; holdvar=prod1m(:,:,i); prodtot1(i,1)=nansum(holdvar(:)); end
for i=1:12; holdvar=prodm(:,:,i); prodtot(i,2)=nansum(holdvar(subtropSPac)); end
for i=1:12; holdvar=prod1m(:,:,i); prodtot1(i,2)=nansum(holdvar(subtropSPac)); end
for i=1:12; holdvar=prodm(:,:,i); prodtot(i,3)=nansum(holdvar(arctic)); end
for i=1:12; holdvar=prod1m(:,:,i); prodtot1(i,3)=nansum(holdvar(arctic)); end
for i=1:12; holdvar=prodm(:,:,i); prodtot(i,4)=nansum(holdvar(osmosis)); end
for i=1:12; holdvar=prod1m(:,:,i); prodtot1(i,4)=nansum(holdvar(osmosis)); end

figure; plot(orderparam12,prodtot1./prodtot,'o'); legend('global','SSP','Arctic','PAP')

figure; plot(orderparam12,(prodtot1-prodtot)./prodtot,'o'); legend('global','SSP','Arctic','PAP')
