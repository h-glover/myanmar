%% High flow - Sept 2017:
% flood is positive, ebb is negative
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_sep17_hc.mat'),hc=aqd;
load('aqd\aqd_sep17_lc.mat'),lc=aqd;clear aqd

% calc instantaneous sed flux in kg/m/s, and fill missing linearly
hc.sf = hc.spd_mean.*hc.ssc1/1000.*hc.depth; % m/s * mg/L / 1000 = (kg/m2/s)*h = kg/s/m
lc.sf = lc.spd_mean.*lc.ssc1/1000.*lc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* h = kg/s/m

hc.sf = fillmissing(hc.sf,'linear');
lc.sf = fillmissing(lc.sf,'linear');

lc.sf(lc.slope<0)=(-1)*lc.sf(lc.slope<0);
hc.sf(hc.slope<0)=(-1)*hc.sf(hc.slope<0);
lc.spd(lc.slope<0,:)=(-1)*lc.spd(lc.slope<0,:);
hc.spd(hc.slope<0,:)=(-1)*hc.spd(hc.slope<0,:);

figure;
subplot(3,2,1)
yyaxis left,pcolor(lc.time,lc.z,lc.spd')
ylim([0 4])
yyaxis right,plot(lc.time,lc.sf,'k-')
ylim([-0.4 0.4])
xlim([lc.time(1) lc.time(end)])
datetick('x','dd','keeplimits')

subplot(3,2,2)
yyaxis left,pcolor(hc.time,hc.z,hc.spd')
ylim([0 6.5])
yyaxis right,plot(hc.time,hc.sf,'k-')
ylim([-2 2])
xlim([hc.time(1) hc.time(end)])
datetick('x','dd','keeplimits')

%%
clear all
load('aqd\aqd_mar18_hc.mat'),hc=aqd;
load('aqd\aqd_mar18_lc.mat'),lc=aqd;clear aqd

% calc instantaneous sed flux in kg/m/s, and fill missing linearly
hc.sf = hc.spd_mean.*hc.ssc1/1000.*hc.depth; % m/s * mg/L / 1000 = (kg/m2/s)*h = kg/s/m
lc.sf = lc.spd_mean.*lc.ssc1/1000.*lc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* h = kg/s/m

hc.sf = fillmissing(hc.sf,'linear');
lc.sf = fillmissing(lc.sf,'linear');

lc.sf(lc.slope<0)=(-1)*lc.sf(lc.slope<0);
hc.sf(hc.slope<0)=(-1)*hc.sf(hc.slope<0);
lc.spd(lc.slope<0,:)=(-1)*lc.spd(lc.slope<0,:);
hc.spd(hc.slope<0,:)=(-1)*hc.spd(hc.slope<0,:);

subplot(3,2,3)
yyaxis left,pcolor(lc.time,lc.z,lc.spd')
ylim([0 4])
yyaxis right,plot(lc.time,lc.sf,'k-')
ylim([-0.4 0.4])
xlim([lc.time(1) lc.time(end)])
datetick('x','dd','keeplimits')

subplot(3,2,4)
yyaxis left,pcolor(hc.time,hc.z,hc.spd')
ylim([0 6.5])
yyaxis right,plot(hc.time,hc.sf,'k-')
ylim([-2 2])
xlim([hc.time(1) hc.time(end)])
datetick('x','dd','keeplimits')

%%
clear all
load('aqd\aqd_sep19_hc.mat'),hc=aqd;
load('aqd\aqd_sep19_lc.mat'),lc=aqd;clear aqd

% calc instantaneous sed flux in kg/m/s, and fill missing linearly
hc.sf = hc.spd_mean.*hc.ssc1/1000.*hc.depth; % m/s * mg/L / 1000 = (kg/m2/s)*h = kg/s/m
lc.sf = lc.spd_mean.*lc.ssc1/1000.*lc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* h = kg/s/m

hc.sf = fillmissing(hc.sf,'linear');
lc.sf = fillmissing(lc.sf,'linear');

lc.sf(lc.slope<0)=(-1)*lc.sf(lc.slope<0);
hc.sf(hc.slope<0)=(-1)*hc.sf(hc.slope<0);

lc.spd(lc.slope<0,:)=(-1)*lc.spd(lc.slope<0,:);
hc.spd(hc.slope<0,:)=(-1)*hc.spd(hc.slope<0,:);

subplot(3,2,5)
yyaxis left,pcolor(lc.time,lc.z,lc.spd')
ylim([0 4])
yyaxis right,plot(lc.time,lc.sf,'k-')
ylim([-0.4 0.4])
xlim([lc.time(1) lc.time(end)])
datetick('x','dd','keeplimits')

subplot(3,2,6)
yyaxis left,pcolor(hc.time,hc.z,hc.spd')
ylim([0 6.5])
yyaxis right,plot(hc.time,hc.sf,'k-')
ylim([-2 2])
xlim([hc.time(1) hc.time(end)])
datetick('x','dd','keeplimits')


%%
for jj=1:6
    subplot(3,2,jj),yyaxis left
    shading flat,colorbar,caxis([-1 1]),colormap(cmocean('balance'))
end