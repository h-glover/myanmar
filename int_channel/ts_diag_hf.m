

%% TS diagrams:
% first look at just internal island ctd casts and aqd casts
% then look at ctds in river

clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('mmi_discharge\ctdprof_13Sep17.mat')
load('aqd\aqd_sep17_hc.mat'),hc=aqd; hc.sal = zeros(size(hc.time));
load('aqd\aqd_sep17_lc.mat'),lc=aqd;clear aqd
load('longterminst\BogaleRiverInstruments.mat')
load('longterminst\BogaleRiverWeather.mat')
load('mmi_discharge\channel_trace.mat')

hc.rho = SW_Density(hc.temp,'C',hc.sal,'ppt',hc.pres./10,'bar');
lc.rho = SW_Density(lc.temp,'C',lc.sal,'ppt',lc.pres./10,'bar');

% generating background density contours
smin = 0;% set min and max values for your plot for sal and temp
smax = 1;
thetamin = 27;
thetamax = 32;
xdim=round((smax-smin)./0.1+1);
ydim=round((thetamax-thetamin)+1);
dens=zeros(ydim,xdim);
thetai=((1:ydim)-1)*1+thetamin;
si=((1:xdim)-1)*0.1+smin;
for j=1:ydim
    for i=1:xdim
        dens(j,i)=SW_Density(thetai(j),'C',si(i),'ppt',0.1,'bar');
    end
end
dens=dens-1000;


%% plotting aqd;
close all
figure;


idx1 = find(hc.time>hc.time(40*60) & hc.time<hc.time(65*60));
idx2 = find(lc.time>hc.time(40*60) & lc.time<hc.time(65*60));
tmin = [(0:1:length(idx1)-1)/60]';

subplot(231)% TS diagrams, colored by slope
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(hc.sal(idx1),hc.temp(idx1),[],hc.slope(idx1),'*')
scatter(lc.sal(idx2),lc.temp(idx2),[],lc.slope(idx2),'o')
colorbar,caxis([-0.8 0.8])
title('aqds with slope'),legend({'dens','HC','LC'})
xlabel('Salinity'),ylabel('Temp (^oC)')
s=gca; s.Colormap=cmocean('balance');

subplot(232)% TS diagrams colored by time:
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(hc.sal(idx1),hc.temp(idx1),[],tmin,'*')
scatter(lc.sal(idx2),lc.temp(idx2),[],tmin,'o')
colorbar,caxis([tmin(1) tmin(end)])
title('aqds with time'),legend({'dens','HC','LC'})
xlabel('Salinity'),ylabel('Temp (^oC)')
s=gca; s.Colormap=cmocean('-thermal');

subplot(233)% stage vel salinity diagrams:
scatter(hc.slope,hc.depth,[],hc.rho-1000,'*'),hold on
scatter(lc.slope,lc.depth,[],lc.rho-1000,'o')
colorbar,caxis([-5 0])
xlabel('slope'),ylabel('depth'),title('aqds with density (kg/m^3)')
s=gca; s.Colormap=cmocean('tempo');

subplot(236)% stage vel ssc diagrams:
scatter(hc.slope,hc.depth,[],hc.ssc2,'*'),hold on
scatter(lc.slope,lc.depth,[],lc.ssc1,'o')
colorbar,caxis([0 200])
xlabel('slope'),ylabel('depth'),title('aqds with ssc (mg/L)')
s=gca; s.Colormap=cmocean('tempo');

subplot(2,3,4:5)
scatter(tmin,hc.depth(idx1),[],hc.slope(idx1),'*'),hold on
scatter(tmin,lc.depth(idx2),[],lc.slope(idx2),'o'),hold on
xlabel(['hour from ',datestr(hc.time(idx1(1)))]),ylabel('depth')
colorbar,caxis([-0.8 0.8])
s=gca; s.Colormap=cmocean('balance');
yyaxis right
plot(tmin(1:10:end),Weather.PAR(Weather.datenum>=hc.time(40*60) & Weather.datenum<=hc.time(65*60)))
ylabel('PAR')
% 1 = lc aqd
% 8 = channel
% 7 = croc
% 6 = turning basin
%% plotting aqd with adcp

idx1 = find(hc.time>hc.time(3850) & hc.time<lc.time(end)); 
% idx1 = 1:length(hc.time);%
idx2 = find(lc.time>hc.time(3850) & lc.time<lc.time(end)); % 
% idx2 = 1:length(lc.time);%


figure;
subplot(3,2,[1,3])
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(hc.sal(idx1),hc.temp(idx1),[],hc.time(idx1),'*')
scatter(lc.sal(idx2),lc.temp(idx2),[],lc.time(idx2),'*')


for jj=1:length(ctd)
%     ctd(jj).Density(ctd(jj).Density<1005)=NaN;
    tvec(jj) = ctd(jj).time(1);
    drho(1,jj)=max(ctd(jj).Density',[],'omitnan');
    drho(2,jj)=min(ctd(jj).Density',[],'omitnan');
    
    subplot(3,2,[1,3])
    scatter(ctd(jj).Salinity,ctd(jj).Temperature,[],ctd(jj).time,'d','filled')
    subplot(3,2,[2,4])
%     scatter(ctd(jj).Density-1000,ctd(jj).Depth,[],ctd(jj).time,'d','filled')
    scatter(ctd(jj).SSCCal,ctd(jj).Depth,[],ctd(jj).time,'d','filled')

    hold on
end
drho(3,:)=drho(1,:)-drho(2,:);
legend('Location','southwest')
axis ij
colorbar,caxis([ctd(1).time(1)-0.04 ctd(end).time(end)])


subplot(3,2,[1,3])
% xlim([0 2]),ylim([27 30])
colorbar,caxis([ctd(1).time(1)-0.04 ctd(end).time(end)])
colormap(cmocean('-thermal'))
title('Aqds with related ctd profiles')
xlabel('Salinity'),ylabel('Temp (^oC)')
legend({'dens','HC','LC','ctd'},'Location','southwest')

subplot(3,2,5:6)
plot(BogaleRiver.datenum,BogaleRiver.FredaDepth,'b'),hold on
plot(hc.time,hc.depth,'k-')
plot(lc.time,lc.depth,'k--')
scatter(tvec,4*ones(size(tvec)),[],tvec,'d','filled')
colorbar,caxis([ctd(1).time(1)-0.04 ctd(end).time(end)])
xlim([ctd(1).time(1)-0.25 ctd(end).time(end)+0.25])
datetick('x','HH','keeplimits')
legend({'Freda','HC','LC','ctd'})


% plot ctd locs along channel:
ctdlat=vertcat(ctd.lat);
ctdlong=vertcat(ctd.long);
figure;
plot(ch.lon,ch.lat,'k-'),hold on
scatter(ctdlong,ctdlat,[],tvec,'d','filled')
caxis([ctd(1).time(1)-0.04 ctd(end).time(end)])
colormap(cmocean('-thermal'))

%%
figure;

idx1 = find(hc.time>hc.time(1) & hc.time<hc.time(360));
idx2 = find(lc.time>hc.time(1) & lc.time<hc.time(360));
tmin = [(0:5:5*357)/60]';

subplot(221)% TS diagrams, colored by slope
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(hc.sal(idx1),hc.temp(idx1),[],hc.slope(idx1),'*')
colorbar,caxis([-0.8 0.8])
title('HC with slope')
xlabel('Salinity'),ylabel('Temp (^oC)')
s=gca; s.Colormap=cmocean('balance');


subplot(222)% TS diagrams, colored by slope
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(lc.sal(idx2),lc.temp(idx2),[],lc.slope(idx2),'o')
colorbar,caxis([-0.8 0.8])
title('LC with slope')
xlabel('Salinity'),ylabel('Temp (^oC)')
s=gca; s.Colormap=cmocean('balance');

subplot(223)% TS diagrams colored by time:
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(hc.sal(idx1),hc.temp(idx1),[],tmin,'*')
colorbar,caxis([tmin(1) tmin(end)])
title('HC with time')
xlabel('Salinity'),ylabel('Temp (^oC)')
s=gca; s.Colormap=cmocean('-thermal');

subplot(224)% TS diagrams colored by time:
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(lc.sal(idx2),lc.temp(idx2),[],tmin,'o')
colorbar,caxis([tmin(1) tmin(end)])
title('LC with time')
xlabel('Salinity'),ylabel('Temp (^oC)')
s=gca; s.Colormap=cmocean('-thermal');

%%

