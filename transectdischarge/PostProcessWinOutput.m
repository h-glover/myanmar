clear all,close all,clc
F=dir('*024*0.mat');
load(F.name)
L=length(AnFLatDeg);

% remove bad values from dir and EN vectors
SerDir10thDeg(abs(SerDir10thDeg)>3600)=NaN;
SerDirDeg=SerDir10thDeg/10;
SerEmmpersec(abs(SerEmmpersec)>15000)=NaN;
SerNmmpersec(abs(SerNmmpersec)>15000)=NaN;
SerMagmmpersec(abs(SerMagmmpersec)>15000)=NaN;

% Correct the EN vectors using nav speed and recalc dir and angle
SerEmmpersec_corr = SerEmmpersec - repmat(AnNVEmmpersec,[1,100]);
SerNmmpersec_corr = SerNmmpersec - repmat(AnNVNmmpersec,[1,100]);
angle_corr = 90-(180/pi)*atan2(SerNmmpersec_corr,SerEmmpersec_corr);
mag_corr = hypot(SerNmmpersec_corr,SerEmmpersec_corr);

%%
figure;
subplot(3,2,1)
plot(AnFLonDeg,AnFLatDeg,'k*'),hold on
plot(AnLLonDeg,AnLLatDeg,'r*')

subplot(3,2,3:4)
plot(AnH100thDeg/100,'r'),ylabel('heading'),hold on
plot(AnNVDir10thDeg/10,'k')
subplot(3,2,5:6)
quiver([1:L]',zeros([L,1]),AnNVEmmpersec,AnNVNmmpersec),ylabel('nav')

%% correct the direction using nav speed, compare to original

figure;
subplot(311)
pcolor(1:L,SerBins(1:25),SerDirDeg(:,1:25)'),shading flat,colorbar
caxis([0 360])

subplot(312)
pcolor(1:L,SerBins(1:25),angle_corr(:,1:25)'),shading flat,colorbar
caxis([0 360])

del = round(SerDirDeg-angle_corr);
del(del>360)=360-del(del>360); del(del==360)=0;
subplot(313)
pcolor(1:L,SerBins(1:25),del(:,1:25)'),shading flat,colorbar
% caxis([0 360])

%% correct the speed using nav speed, compare to original

figure;
subplot(311)
pcolor(1:L,SerBins(1:25),[SerMagmmpersec(:,1:25)/1000]')
shading flat,colorbar, caxis([0 3])

subplot(312)
pcolor(1:L,SerBins(1:25),[mag_corr(:,1:25)/1000]')
shading flat,colorbar,caxis([0 3])

del = SerMagmmpersec-mag_corr;
subplot(313)
pcolor(1:L,SerBins(1:25),[del(:,1:25)/1000]'),shading flat,colorbar
caxis([-1 1])
%%
% angle = 90-(180/pi)*atan2(adcp.north_vel,adcp.east_vel);
% u = EW
% v = NS
% angle_corr=angle_corr(20:end,1:36);
% mag_corr=mag_corr(20:end,1:36);

%% look at the E and N corrected values (should be flooding (~290))

figure;
subplot(211)
pcolor(1:L,SerBins(1:25),[SerEmmpersec_corr(:,1:25)/1000]')
shading flat,colorbar,caxis([-5 5])
title('East')

subplot(212)
pcolor(1:L,SerBins(1:25),[SerNmmpersec_corr(:,1:25)/1000]')
shading flat,colorbar,caxis([-5 5])
title('North')

figure;
subplot(211)
pcolor(1:L,SerBins(1:25),[SerEmmpersec(:,1:25)/1000]')
shading flat,colorbar,caxis([-5 5])
title('East,uncorr')

subplot(212)
pcolor(1:L,SerBins(1:25),[SerNmmpersec(:,1:25)/1000]')
shading flat,colorbar,caxis([-5 5])
title('North,uncorr')
