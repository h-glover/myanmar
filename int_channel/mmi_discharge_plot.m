%% plort
clear all,close all,clc
cd C:\GLOVER\output\myanmar

load('mmi_discharge\br_int_march2018.mat')
load('longterminst\BogaleRiverInstruments.mat'),br=BogaleRiver; clear BogaleRiver

%interp to a uniform depth interval
intrp_depth =0:0.05:15; %max adcp depth is 15m

for kk=[3]
figure;
subplot(3,3,1:2)
pcolor(adcp(kk).elapdist,adcp(kk).z,adcp(kk).spd)
shading flat,colorbar,ylim([0 13]),caxis([0 50]),title('speed'),axis ij
subplot(3,3,4:5)
pcolor(adcp(kk).elapdist,intrp_depth,adcp(kk).ssc)
shading flat,colorbar,ylim([0 13]),caxis([0 1000]),title('SSC'),axis ij
subplot(3,3,7:8)
pcolor(adcp(kk).elapdist,intrp_depth,adcp(kk).sal)
shading flat,colorbar,ylim([0 13]),caxis([15 20]),title('Salinity'),axis ij
subplot(3,3,[3,6])
% quiver(adcp(kk).lon,adcp(kk).lat,adcp(kk).east(2,:),adcp(kk).north(2,:))
scatter(adcp(kk).lon,adcp(kk).lat,[],adcp(kk).elapdist),colorbar
subplot(3,3,9)
yyaxis left
plot(adcp(kk).time,adcp(kk).elapdist,'k-')%,ylim([-5000 10000])
yyaxis right
plot(br.datenum,br.ut_FredaDepth,'b-')
xlim([adcp(1).time(1) adcp(3).time(end)])
colormap(cmocean('amp'))
figure;
ee = nansum(adcp(kk).sed_flux_east);%nanmean(adcp(kk).sed_flux_east);
ee = movmean(ee,5);
nn = nansum(adcp(kk).sed_flux_north);%nanmean(adcp(kk).sed_flux_north);
nn = movmean(nn,5);
subplot(121)
quiver(adcp(kk).lon,adcp(kk).lat,ee,nn,'AutoScaleFactor',10,'ShowArrowHead','off')
title('sediment flux')

% ee = nansum(adcp(kk).sal_flux_east);
% ee = movmean(ee,5);
% nn = nansum(adcp(kk).sal_flux_north);
% nn = movmean(nn,5);
% subplot(132)
% quiver(adcp(kk).lon,adcp(kk).lat,ee,nn,'AutoScaleFactor',5,'ShowArrowHead','off')
% title('salt flux')

ee = nansum(adcp(kk).east);
ee = movmean(ee,5);
nn = nansum(adcp(kk).north);
nn = movmean(nn,5);
subplot(122)
quiver(adcp(kk).lon,adcp(kk).lat,ee,nn,'AutoScaleFactor',5,'ShowArrowHead','off')
title('velocity')
end



%% fig for S2S talk:

clear all,close all,clc
cd C:\GLOVER\output\myanmar

load('mmi_discharge\br_int_march2018.mat')
load('longterminst\BogaleRiverInstruments.mat'),br=BogaleRiver; clear BogaleRiver

kk=3;
%interp to a uniform depth interval
intrp_depth =0:0.05:15; %max adcp depth is 15m
adcp(kk).elapdist = adcp(kk).elapdist/1000;

figure;
subplot(3,1,1)
pcolor(adcp(kk).elapdist,adcp(kk).z,adcp(kk).spd)
caxis([0 50])
ax=gca;ax.Colormap = cmocean('amp');
subplot(3,1,2)
pcolor(adcp(kk).elapdist,intrp_depth,adcp(kk).ssc)
caxis([0 1000]),ax=gca;ax.Colormap = cmocean('turbid');
subplot(3,1,3)
pcolor(adcp(kk).elapdist,intrp_depth,adcp(kk).sal)
caxis([15 20]),ax=gca;ax.Colormap = cmocean('haline');



for jj=1:3
    subplot(3,1,jj)
    shading flat
    colorbar
    ylim([0 13])
    axis ij
    ylabel('Depth (m)'),xlabel('Dist. from southern entrance (km)')
end
