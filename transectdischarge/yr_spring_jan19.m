%% Yangon Spring Low Flow Data:

clear all%,close all,clc
cd C:\GLOVER\data\myanmar\transecdischarge\YR_Spring_Jan19
% 
% % load water level data:
% load('C:\GLOVER\output\myanmar\longterminst\YangonRiverInstruments.mat')

load('ENX29')
enx.bindepth=RDIBin1Mid:RDIBinSize:(RDIBinSize*100+RDIBin1Mid/2);
L=length(SerHour);
enx.tvec=datenum(SerYear,SerMon,SerDay,SerHour,SerMin,SerSec);
[enx.boatspd,enx.boatdir]=uv2sd(...
    AnNVEmmpersec,AnNVNmmpersec,zeros(length(AnNVEmmpersec),1));
enx.v1=SerEmmpersec/1000; 
enx.v2=SerNmmpersec/1000;


load('ENS29')
ens.v1=SerEmmpersec/1000; 
ens.v2=SerNmmpersec/1000; 

figure;scatter(ens.v1(:,2),enx.v1(:,2))
figure;scatter(ens.v2(:,2),enx.v2(:,2))

%%
figure;
quiver(AnFLonDeg,AnFLatDeg,AnNVEmmpersec,AnNVNmmpersec,'k'),hold on
quiver(AnFLonDeg,AnFLatDeg,SerEmmpersec(:,2),SerNmmpersec(:,2),'r')
%%
% remove boat velocity from water velocity
SerEmmpersec=SerEmmpersec-AnNVEmmpersec;
SerNmmpersec=SerNmmpersec-AnNVEmmpersec;

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


% clearvars -except enx
% load('ENS29.mat')

% 1.4797   -1.4338   -0.0797    0.0454
% 0.0486   -0.0770   -1.4308    1.4692
% 0.2588    0.2740    0.2727    0.259  
% 0.9987    1.0575   -1.0537   -1.0038


% % enu conversion for within burst velocity
% % use burst avg head, pitch, roll and compute for each burst
% for jj=1:L
%     [vec.v1b(jj,:),vec.v2b(jj,:),vec.v3b(jj,:)]=xyz2enu(...
%         vec.v1b(jj,:),vec.v2b(jj,:),vec.v3b(jj,:),...
%         vec.head(jj),vec.pitch(jj),vec.roll(jj),vec.meta.transmat/4096);
% end
