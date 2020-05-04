% down river is 105 (almost directly east)
% up river is 285 (almost directly west)
% btwn 15-195=down river


clear all,close all,clc

v=VideoWriter('YR_Jan19_Spring_QC');
v.FrameRate=1;
open(v)

F=dir('ADCP*.mat');
C1=flipud([ones([30,1]),linspace(1,0,30)',linspace(1,0,30)']);
C2=[linspace(1,0,30)',linspace(1,0,30)',ones([30,1])];
C=[C1(1:end-3,:);C2(4:end,:)];

for jj=1:length(F)
    load(F(jj).name)
    blankdist=RDIBin1Mid-RDIBinSize;
    depth=RDIBin1Mid:RDIBinSize:(RDIBinSize*length(SerBins)+blankdist);
    
    SerMagmmpersec=SerMagmmpersec'/1000;
    SerMagmmpersec(SerMagmmpersec<0 | SerMagmmpersec>4)=NaN;
    SerDir10thDeg=SerDir10thDeg'/10;
    SerDir10thDeg(SerDir10thDeg<0)=NaN;
    
    [m,n]=size(SerDir10thDeg);
    dwnrvr=(-1)*ones(m,n); %make river dir -1 for up river and 1 for down river
    dwnrvr(SerDir10thDeg>15 & SerDir10thDeg<195)=1;
    dwnrvr(isnan(SerDir10thDeg))=NaN;
    
    SerEmmpersec=SerEmmpersec'/1000;
    SerEmmpersec(abs(SerEmmpersec)>4)=NaN;
    
    t(jj)=datenum([RDIEnsDate,RDIEnsTime(1:end-6)],'yy/mm/ddHH:MM');
    spd(jj)=nanmean(nanmean(SerMagmmpersec,1));
    evel(jj)=nanmean(nanmean(SerEmmpersec,1));
    
    
    fig=figure('visible','off');
    subplot(311)
    pcolor(1:SerEnsembles(end),depth,SerMagmmpersec)
    shading flat
    ylim([0 22]),axis ij,ylabel('Speed (m/s)')
    colorbar,caxis([0 4])
    title(RDIEnsTime)
    
    subplot(312)
    pcolor(1:SerEnsembles(end),depth,SerDir10thDeg)
    shading flat
    ylim([0 22]),axis ij,ylabel('Direction')
    colorbar,caxis([0 360])
    
    subplot(313)
    pcolor(1:SerEnsembles(end),depth,SerEmmpersec)
    shading flat
    ylim([0 22]),axis ij,ylabel('East Velocity (m/s)')
    colorbar,caxis([-4 4]),colormap(C)
    
    fg(jj)=getframe(fig,[0 0 560 420]);
    writeVideo(v,fg(jj))
    
    maxeast(jj,1)=max(max(abs(SerEmmpersec)));
    maxspd(jj,1)=max(max(abs(SerMagmmpersec)));
    
end
close(v)
close all

load('YangonWL_Prediction.mat')
p_time2=days([min(p_time):minutes(1):max(p_time)]);
y_out2=interp1(p_time,YOUT2,p_time2);

[~,iwl,iA]=intersect(p_time2,t);
wl=y_out2(iwl)';

figure;plot(p_time,YOUT2)
xlim([t(1) t(end)]),datetick('x','HH','keeplimits')
figure;
plot(maxeast,'r'),hold on,plot(maxspd,'k'),legend('east','spd')


%% check lat/long boat position data:

clear all,close all,clc

F=dir('ADCP*.mat');
leftbank=[16.766193,96.17034];
rightbank=[16.75907,96.167959];
for jj=1:length(F)
    load(F(jj).name)
    t(jj)=datenum([RDIEnsDate,RDIEnsTime(1:end-6)],'yy/mm/ddHH:MM');
    
    figure;
    subplot(311)
    quiver(AnFLonDeg(20:end),AnFLatDeg(20:end),AnNVEmmpersec(20:end)./100,AnNVNmmpersec(20:end)./100)
    hold on
    plot(leftbank(2),leftbank(1),'r*'),plot(rightbank(2),rightbank(1),'g*')
    ylabel('lat/lon optns')
    title(RDIEnsTime)
    
    subplot(312)
    plot(abs(AnH100thDeg)/100),hold on,plot(AnNVDir10thDeg/10)
    legend('Head','NavDir')
    ylabel('Heading/NavDir')
    
    subplot(313)
    plot(abs(AnH100thDeg)/100-AnNVDir10thDeg/10)
end

%% test the 30th transect to see if navigation data makes sense
clear all,clc
load('test30.mat')
leftbank=[16.766193,96.17034];
rightbank=[16.75907,96.167959];

figure;
subplot(311)
quiver(AnFLonDeg(20:end),AnFLatDeg(20:end),AnNVEmmpersec(20:end)./100,AnNVNmmpersec(20:end)./100)
hold on
plot(leftbank(2),leftbank(1),'r*'),plot(rightbank(2),rightbank(1),'g*')
ylabel('lat/lon optns')
title(RDIEnsTime)

subplot(312)
plot(abs(AnH100thDeg(20:end))/100),hold on,plot(AnNVDir10thDeg(20:end)/10)
legend('Head','NavDir')
ylabel('Heading/NavDir')

subplot(313)
plot(abs(AnH100thDeg(20:end))/100-AnNVDir10thDeg(20:end)/10)
xlabel('HDG v HDT')


%% test the 30th transect, referenced to VMDAS NAV in WinADCP
clear all,close all,clc
C1=flipud([ones([30,1]),linspace(1,0,30)',linspace(1,0,30)']);
C2=[linspace(1,0,30)',linspace(1,0,30)',ones([30,1])];
C=[C1(1:end-3,:);C2(4:end,:)];

load('test24.mat')
blankdist=RDIBin1Mid-RDIBinSize;
depth=RDIBin1Mid:RDIBinSize:(RDIBinSize*length(SerBins)+blankdist);

SerMagmmpersec=SerMagmmpersec'/1000;
SerMagmmpersec(SerMagmmpersec<0 | SerMagmmpersec>4)=NaN;
SerDir10thDeg=SerDir10thDeg'/10;
SerDir10thDeg(SerDir10thDeg<0)=NaN;

[m,n]=size(SerDir10thDeg);
dwnrvr=(-1)*ones(m,n); %make river dir -1 for up river and 1 for down river
dwnrvr(SerDir10thDeg>15 & SerDir10thDeg<195)=1;
dwnrvr(isnan(SerDir10thDeg))=NaN;

SerEmmpersec=SerEmmpersec'/1000;
SerEmmpersec(abs(SerEmmpersec)>4)=NaN;

meanspd=nanmean(nanmean(SerMagmmpersec,1))
meanevel=nanmean(nanmean(SerEmmpersec,1))
maxspd=max(max(abs(SerMagmmpersec)))
maxeast=max(max(abs(SerEmmpersec)))

figure;
subplot(311)
pcolor(1:SerEnsembles(end),depth,SerMagmmpersec)
shading flat
ylim([0 22]),axis ij,ylabel('Speed (m/s)')
colorbar,caxis([0 4])
title(RDIEnsTime)

subplot(312)
pcolor(1:SerEnsembles(end),depth,SerDir10thDeg)
%pcolor(1:SerEnsembles(end),depth,dwnrvr)
shading flat
ylim([0 22]),axis ij,ylabel('Direction')
colorbar,caxis([-1 1]),caxis([0 360])

subplot(313)
pcolor(1:SerEnsembles(end),depth,SerEmmpersec)
shading flat
ylim([0 22]),axis ij,ylabel('East Velocity (m/s)')
colorbar,caxis([-4 4]),%colormap(C)
