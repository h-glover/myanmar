%% what does yangon neap low flow data look like at similar point in tide?

clear all,close all,clc
cd C:\GLOVER\data\myanmar\transecdischarge\YR_Neap_Jan19
% % load water level data:
load('C:\GLOVER\output\myanmar\longterminst\YangonRiverInstruments.mat')
figure;
subplot(211)
plot(YangonRiver.datenum,YangonRiver.ut_Depth),hold on
idx=find(YangonRiver.datenum==datenum('21-Jan-2019 17:20'));
plot(YangonRiver.datenum(idx),YangonRiver.ut_Depth(idx),'ro')
xlim([datenum('21-Jan-2019') datenum('22-Jan-2019')])
datetick('x','HH','keeplimits')
subplot(212)
plot(YangonRiver.datenum,YangonRiver.ut_Depth),hold on
xlim([datenum('14-Jan-2019') datenum('15-Jan-2019')])
datetick('x','HH','keeplimits')

load('YR_Jan19_Neap.mat')
for jj=1:length(adcp)
t(jj)=adcp(jj).time(1);
end
plot(t,zeros(1,length(t)),'ro')

% first/2nd transects are most comparable to 23 (towards end of flood)

jj=7;
% figure;
% quiver(adcp(jj).lon,adcp(jj).lat,adcp(jj).east(2,:),adcp(jj).north(2,:))

[spd,dir]=uv2sd(adcp(jj).east,adcp(jj).north,adcp(jj).up);

figure;
subplot(211),pcolor(adcp(jj).time,adcp(jj).z(:,jj),spd)
shading flat,colorbar,caxis([0 150])
ylim([0 20]),axis ij
subplot(212),pcolor(adcp(jj).time,adcp(jj).z(:,jj),dir)
shading flat,colorbar,caxis([0 360])
ylim([0 20]),axis ij
yyaxis right,plot(adcp(jj).time,nanmean(dir),'k','LineWidth',2)
ylim([ 0 360])

%% Yangon Spring Low Flow Data:

clear all%,close all,clc
cd C:\GLOVER\data\myanmar\transecdischarge\YR_Spring_Jan19
% 

load('ENX24')
enx.bindepth=RDIBin1Mid:RDIBinSize:(RDIBinSize*100+RDIBin1Mid/2);
L=length(SerHour);
enx.tvec=datenum(SerYear,SerMon,SerDay,SerHour,SerMin,SerSec);
[enx.boatspd,enx.boatdir]=uv2sd(...
    AnNVEmmpersec,AnNVNmmpersec,zeros(length(AnNVEmmpersec),1));

% % remove boat velocity from water velocity
% SerEmmpersec=SerEmmpersec-AnNVEmmpersec;
% SerNmmpersec=SerNmmpersec-AnNVEmmpersec;

% convert to m/s and remove high vals
enx.v1=SerEmmpersec/1000; % beam 1
enx.v2=SerNmmpersec/1000; % beam 2
enx.v3=SerVmmpersec/1000; % beam 3
enx.v1(abs(enx.v1)>4)=NaN;
enx.v2(abs(enx.v2)>4)=NaN;
enx.v3(abs(enx.v3)>4)=NaN;
enx.v1(:,enx.bindepth>20)=NaN;
enx.v2(:,enx.bindepth>20)=NaN;
enx.v3(:,enx.bindepth>20)=NaN;

[enx.spd,enx.dir]=uv2sd(enx.v1,enx.v2,enx.v3);

% high flow speed rarely exceeds 2.5 m/s
enx.spd(enx.spd>2.5)=NaN;

figure;
subplot(211),pcolor(enx.tvec,enx.bindepth,enx.spd')
shading flat,colorbar,caxis([0 1.5])
ylim([0 20]),axis ij
subplot(212),pcolor(enx.tvec,enx.bindepth,enx.dir')
shading flat,colorbar,caxis([0 360])
ylim([0 20]),axis ij
yyaxis right,plot(enx.tvec,nanmean(enx.dir,2),'k','LineWidth',2)
ylim([ 0 360])

%% compare with return transect:
clear all%,close all,clc
cd C:\GLOVER\data\myanmar\transecdischarge\YR_Spring_Jan19
% 

load('ENX24')
enx.bindepth=RDIBin1Mid:RDIBinSize:(RDIBinSize*100+RDIBin1Mid/2);
L=length(SerHour);
enx.tvec=datenum(SerYear,SerMon,SerDay,SerHour,SerMin,SerSec);
[enx.boatspd,enx.boatdir]=uv2sd(...
    AnNVEmmpersec,AnNVNmmpersec,zeros(length(AnNVEmmpersec),1));

% % remove boat velocity from water velocity
% SerEmmpersec=SerEmmpersec-AnNVEmmpersec;
% SerNmmpersec=SerNmmpersec-AnNVEmmpersec;

% convert to m/s and remove high vals
enx.v1=SerEmmpersec/1000; % beam 1
enx.v2=SerNmmpersec/1000; % beam 2
enx.v3=SerVmmpersec/1000; % beam 3
enx.v1(abs(enx.v1)>4)=NaN;
enx.v2(abs(enx.v2)>4)=NaN;
enx.v3(abs(enx.v3)>4)=NaN;
enx.v1(:,enx.bindepth>20)=NaN;
enx.v2(:,enx.bindepth>20)=NaN;
enx.v3(:,enx.bindepth>20)=NaN;

[enx.spd,enx.dir]=uv2sd(enx.v1,enx.v2,enx.v3);

% high flow speed rarely exceeds 2.5 m/s
enx.spd(enx.spd>2.5)=NaN;

figure;
subplot(211),pcolor(enx.tvec,enx.bindepth,enx.spd')
shading flat,colorbar,caxis([0 1.5])
ylim([0 20]),axis ij
subplot(212),pcolor(enx.tvec,enx.bindepth,enx.dir')
shading flat,colorbar,caxis([0 360])
ylim([0 20]),axis ij
yyaxis right,plot(enx.tvec,nanmean(enx.dir,2),'k','LineWidth',2)
ylim([ 0 360])





%%

% %%
% figure(11)
% quiver(adcp(jj).lon,adcp(jj).lat,adcp(jj).east(2,:),adcp(jj).north(2,:)),hold on
% figure(11)
% quiver(AnFLonDeg,AnFLatDeg,enx.v1(:,2),enx.v2(:,2),'r--')
% figure(11)
% quiver(AnFLonDeg,AnFLatDeg,enx.v1(:,2),enx.v2(:,2),'b--')
% legend({'adcp','23','24'})

% rb = [96.166544,16.759529];
% lb = [96.168025,16.766170];
%
% figure;
% subplot(311),pcolor(enx.tvec,enx.bindepth,enx.v1'),shading flat,colorbar,ylim([0 20])
% subplot(312),pcolor(enx.tvec,enx.bindepth,enx.v2'),shading flat,colorbar,ylim([0 20])
% subplot(313),pcolor(enx.tvec,enx.bindepth,enx.v3'),shading flat,colorbar,ylim([0 20])
% 
% figure;
% scatter(AnFLonDeg,AnFLatDeg,[],boatdir),colorbar,hold on,
% quiver(AnFLonDeg,AnFLatDeg,AnNVEmmpersec,AnNVNmmpersec,'k')
% plot(lb(1),lb(2),'ro'),plot(rb(1),rb(2),'go')



% % enu conversion for within burst velocity
% % use burst avg head, pitch, roll and compute for each burst
% for jj=1:L
%     [vec.v1b(jj,:),vec.v2b(jj,:),vec.v3b(jj,:)]=xyz2enu(...
%         vec.v1b(jj,:),vec.v2b(jj,:),vec.v3b(jj,:),...
%         vec.head(jj),vec.pitch(jj),vec.roll(jj),vec.meta.transmat/4096);
% end


% load('ENS29')
% ens.v1=SerEmmpersec/1000; 
% ens.v2=SerNmmpersec/1000; 
% ens.boat1=AnNVEmmpersec/1000;
% ens.boat2=AnNVEmmpersec/1000;
% ens.sub1=ens.v1-ens.boat1;
% ens.sub2=ens.v2-ens.boat2;
% 
% figure;
% quiver(AnFLonDeg,AnFLatDeg,ens.boat1,ens.boat2,'k'),hold on
% quiver(AnFLonDeg,AnFLatDeg,ens.v1(:,2),ens.v2(:,2),'r')
% quiver(AnFLonDeg,AnFLatDeg,ens.sub1(:,2),ens.sub2(:,2),'b')
% quiver(AnFLonDeg,AnFLatDeg,enx.v1(:,2),enx.v2(:,2),'r--')
% quiver(AnFLonDeg,AnFLatDeg,enx.sub1(:,2),enx.sub2(:,2),'b--')
% 
% %
% figure;
% quiver(AnFLonDeg,AnFLatDeg,AnNVEmmpersec,AnNVNmmpersec,'k'),hold on
% quiver(AnFLonDeg,AnFLatDeg,SerEmmpersec(:,2),SerNmmpersec(:,2),'r')