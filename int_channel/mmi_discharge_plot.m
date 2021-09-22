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

%%
%% 
clear all%,close all,clc

cd C:\GLOVER\output\myanmar
load('mmi_discharge\br_int_march2018.mat')
load('longterminst\BogaleRiverInstruments.mat')
br=BogaleRiver;

for jj=1:3
adcp(jj).east_mean = nanmean(adcp(jj).east);
adcp(jj).north_mean = nanmean(adcp(jj).north);

adcp(jj).heading(adcp(jj).heading<0)=NaN;
adcp(jj).dir(adcp(jj).dir<0)=NaN;
adcp(jj).spd(adcp(jj).spd>200)=NaN;

adcp(jj).spd_mean = nanmean(adcp(jj).spd);
adcp(jj).spd_mean = movmean(adcp(jj).spd_mean,7);
adcp(jj).dir_mean = nanmean(adcp(jj).dir);
adcp(jj).dir_mean = movmean(adcp(jj).dir_mean,7);


figure(1);
subplot(211)
yyaxis left
plot(adcp(jj).time,adcp(jj).spd_mean,'k-'),hold on
yyaxis right
plot(br.datenum,br.ut_FredaDepth,'b-')
xlim([adcp(1).time(1) adcp(3).time(end)])
datetick('x','HH:MM','keeplimits')
% title(datestr(floor(adcp(jj).time(1))))
hold on
subplot(212)
yyaxis left
plot(adcp(jj).time,adcp(jj).dir_mean,'k-'),hold on
plot(adcp(jj).time,adcp(jj).heading,'r-'),hold on
yyaxis right
plot(br.datenum,br.ut_FredaDepth,'b-')
xlim([adcp(1).time(1) adcp(3).time(end)])
datetick('x','HH:MM','keeplimits')
hold on

% figure(2);
% scatter(adcp(jj).dir_mean,adcp(jj).heading,[],adcp(jj).time),hold on
% colorbar,refline(1,0),refline(-1,350)

adcp(jj).flowdiff = abs(adcp(jj).dir_mean - adcp(jj).heading);
figure(2);
scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).flowdiff,'.')
colorbar,caxis([0 180]),hold on,title('water spd')
colormap(cmocean('phase')),title('water dir - head')

figure;%(2)
scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).spd_mean)
colorbar,caxis([0 50]),hold on,title('water spd')
figure;%(3)
subplot(121),scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).dir_mean)
colorbar,caxis([0 360]),hold on,title('water dir')
subplot(122),scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).heading)
colorbar,caxis([0 360]), hold on
colormap(cmocean('phase')),title('head')
end


%%
clear all,close all,clc

cd C:\GLOVER\output\myanmar\mmi_discharge
load('br_int_march2018')

for jj=1:3
adcp(jj).dir(adcp(jj).dir<0)=NaN;
adcp(jj).spd(adcp(jj).spd>200)=NaN;

figure(10);
quiver(adcp(jj).lon,adcp(jj).lat,adcp(jj).east(2,:),adcp(jj).north(2,:))
hold on,

figure(20);
scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).spd(2,:))
hold on,colorbar,caxis([0 50]),title('water speed')

figure(30);
scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).depth(1,:)),colorbar
hold on,colorbar,title('depth')
end


%% fig for paper: LF velocity, turb, sal, Ri

clear all%,close all,clc
cd C:\GLOVER\output\myanmar

load('mmi_discharge\br_int_march2018.mat')

adcp(3).Ri(:,3) = adcp(3).Ri(:,4);
adcp(3).Ri(1:4,:) = fillmissing(adcp(3).Ri(1:4,:),'nearest',2);
% adcp(3).Ri = movmean(adcp(3).Ri,3,2);
% adcp(3).Ri = movmean(adcp(3).Ri,3);
for jj=1:length(adcp(3).time)
    adcp(3).Ri(adcp(3).z(:,1)>adcp(3).depth(jj),jj) = NaN;
end   
% K = (1/9)*ones(3);
% adcp(3).Ri = conv2(adcp(3).Ri,K,'same');

%interp to a uniform depth interval
intrp_depth = 0:0.05:15; %max adcp depth is 15m
adcp(3).elapdist = adcp(3).elapdist/1000;

figure;
subplot(4,2,1)
pcolor(adcp(3).elapdist,adcp(3).zComplete,adcp(3).spdComplete)
caxis([0 50])
ax=gca;ax.Colormap = cmocean('amp');title('velocity (cm/s)')
subplot(4,2,3)
pcolor(adcp(3).elapdist,adcp(3).z(:,1),adcp(3).Ri)
% pcolor(adcp(3).elapdist,adcp(3).zComplete(:,1),adcp(3).Ri)
caxis([0 1]),ax = gca; ax.Colormap = cmocean('balance','pivot',0.25);title('Ri')
subplot(4,2,5)
pcolor(adcp(3).elapdist,intrp_depth,adcp(3).sal)
caxis([15 20]),ax=gca;ax.Colormap = cmocean('haline');title('Salinity')
subplot(4,2,7)
pcolor(adcp(3).elapdist,intrp_depth,adcp(3).ssc)
caxis([0 1000]),ax=gca;ax.Colormap = cmocean('turbid');title('SSC (mg/L)')

for jj=[1,3,5,7]
    subplot(4,2,jj)
    shading interp
    colorbar
    ylim([0 13]),xlim([0 14.5])
    axis ij
    ylabel('Depth (m)')
end
xlabel('Dist. from southern entrance (km)')


load('mmi_discharge\ctdprof_13Sep17.mat')

ctd_locs = vertcat(ctd.dist);
ctd_locs = round(ctd_locs);
adcp_dist = adcp(3).elapdist*1000;
[adcp_dist,idx]=unique(adcp_dist);
adcp_depth = adcp(3).depth(idx);
dist = 1:max(ctd_locs);
adcp_depth = interp1(adcp_dist,adcp_depth,dist);

ssc = NaN(length(intrp_depth),length(dist));% ssc = NaN(length(intrp_depth),length(ctd));
sal = ssc;
depth =NaN(1,length(dist));

for jj=1:length(ctd)
    ssc(:,dist==ctd_locs(jj)) = interp1(ctd(jj).Depth,ctd(jj).SSCCal,intrp_depth)';
    sal(:,dist==ctd_locs(jj)) = interp1(ctd(jj).Depth,ctd(jj).Salinity,intrp_depth)';
%     depth(1,dist==ctd_locs(jj)) = max(ctd(jj).Depth);
%     ssc(:,jj) = interp1(ctd(jj).Depth,ctd(jj).SSCCal,intrp_depth)';
%     sal(:,jj) = interp1(ctd(jj).Depth,ctd(jj).Salinity,intrp_depth)';
end
% depth=fillmissing(depth,'linear');
ssc = fillmissing(ssc,'nearest',1);
sal = fillmissing(sal,'nearest',1);
ssc = fillmissing(ssc,'linear',2);
sal = fillmissing(sal,'linear',2);
% nan values below the riverbed
for jj=1:length(depth)
    ssc(intrp_depth>adcp_depth(jj),jj) = NaN;
    sal(intrp_depth>adcp_depth(jj),jj) = NaN;
end

subplot(4,2,6)
% pcolor(ctd_locs/1000,intrp_depth,sal)
pcolor(dist/1000,intrp_depth,sal)
caxis([0 1]),ax=gca;ax.Colormap = cmocean('haline');title('Salinity')
subplot(4,2,8)
pcolor(dist/1000,intrp_depth,ssc)
caxis([0 1000]),ax=gca;ax.Colormap = cmocean('turbid');title('SSC (mg/L)')

for jj=[6,8]
    
    
    subplot(4,2,jj)
    shading interp
    colorbar
    ylim([0 13]),xlim([0 14.5])
    axis ij
end
xlabel('Dist. from southern entrance (km)')
clear all

% for jj=[1,3,5:8]
%     subplot(4,2,jj)
%     cla
% end