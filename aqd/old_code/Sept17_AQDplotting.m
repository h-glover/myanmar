% aqd plotting
clear all,close all,clc
load('MMsm_Sept17_aqd')
smch=aqd;
load('MMLG_Sept17_aqd')
lgch=aqd; clear aqd

% take value every 5 min for sm channel record
smch.time=smch.time(1:5:end);
smch.pres=smch.pres(1:5:end);
smch.ext1=smch.ext1(1:5:end);
smch.sal=smch.sal(1:5:end);
lgch.time=lgch.time(1:5:end);
lgch.pres=lgch.pres(1:5:end);
lgch.ext1=lgch.ext1(1:5:end);
lgch.rbr.salinity=lgch.rbr.salinity(1:5:end);

% calculate SSC using Mekong conversion
lgch.ssc1=lgch.ext1.*0.0688 + 50;
lgch.rbr.ssc=lgch.rbr.turb*0.0015*1000;
smch.ssc1=smch.ext1.*0.0688 + 50;

% calculate slope of pressure in m/hr
lgch.slope=lgch.pres(2:end)-lgch.pres(1:end-1);
lgch.slope=[lgch.slope(1);lgch.slope]*12;
smch.slope=smch.pres(2:end)-smch.pres(1:end-1);
smch.slope=[smch.slope(1);smch.slope]*12;

% find values when instruments are in water
smidx=find(smch.pres>0.2);
lgidx=find(lgch.pres>1);
lgidx=lgidx(355:end);
%% Plot the time series of pres, temp, sal, ssc
figure;
subplot(411)
plot(smch.time(smidx),smch.pres(smidx),'b'),hold on
plot(lgch.time(lgidx),lgch.pres(lgidx),'k')
ylim([0 7]),ylabel('Depth (m)')
xlim([lgch.time(lgidx(1)) smch.time(smidx(end))])
datetick('x','dd','keeplimits')
legend('Dead-end Channel','Main Channel')
%
subplot(412)
plot(smch.time(smidx),smch.temp(smidx),'b'),hold on
plot(lgch.time(lgidx),lgch.temp(lgidx),'k')
xlim([lgch.time(lgidx(1)) smch.time(smidx(end))])
ylim([27 32]),ylabel('Temp (C)')
datetick('x','dd','keeplimits')
%
subplot(413)
plot(smch.time(smidx),smch.sal(smidx),'b'),hold on
plot(lgch.time(lgidx),lgch.rbr.salinity(lgidx),'k')
xlim([lgch.time(lgidx(1)) smch.time(smidx(end))])
ylim([0 1.5]),ylabel('Salinity')
datetick('x','dd','keeplimits')

subplot(414)
plot(smch.time(smidx),smch.ssc1(smidx),'b'),hold on
plot(lgch.time(lgidx),lgch.ssc1(lgidx),'k')
xlim([lgch.time(lgidx(1)) smch.time(smidx(end))])
datetick('x','dd','keeplimits')
ylim([0 4000]),ylabel('SSC (mg/L) at 20 cmab')
xlabel(['Day in ',datestr(lgch.time(1),'mmm ''yy')])



%% plot tidal slope v salinity v ssc (color)
figure;
subplot(2,2,1)
scatter(lgch.slope(lgidx),lgch.rbr.salinity(lgidx),[],lgch.ssc1(lgidx),'.')
caxis([0 700])
colorbar
C=[0.8*ones([length(lgidx),1]),linspace(0.8,0,length(lgidx))',linspace(0.8,0,length(lgidx))'];
colormap(C)
axis([-1 1 0 7])
title('High Flow: Large Channel')
ylabel('Salinity')
subplot(2,2,2)
scatter(smch.slope(smidx),smch.sal(smidx),[],smch.ssc1(smidx),'.')
colorbar,caxis([0 700])
axis([-1 1 0 7])
colormap(C)
title('High Flow: Dead-end Channel')

%% salinity v ssc with slope as color
C1=flipud([ones([30,1]),linspace(1,0,30)',linspace(1,0,30)']);
C2=[linspace(1,0,30)',linspace(1,0,30)',ones([30,1])];
C=[C1(1:end-4,:);C2(5:end,:)];

figure;
subplot(2,2,1)
scatter(lgch.rbr.salinity(lgidx),lgch.ssc1(lgidx),[],lgch.slope(lgidx),'.')
caxis([-1 1]),colorbar,colormap(C)
%axis([0 7 0 2000])
title('High Flow: Large Channel')
ylabel('SSC')
subplot(2,2,2)
scatter(smch.sal(smidx),smch.ssc1(smidx),[],smch.slope(smidx),'.')
caxis([-1 1]),colorbar
%axis([0 7 0 2000])
title('High Flow: Dead-end Channel')


%% Just plot salinity v ssc OR slope v ssc
figure;
subplot(2,2,1)
scatter(lgch.rbr.salinity(lgidx),lgch.ssc1(lgidx),'.')
axis([0 7 0 2000])
title('High Flow: Large Channel')
ylabel('SSC')
subplot(2,2,2)
scatter(smch.sal(smidx),smch.ssc1(smidx),'.')
axis([0 7 0 2000])
title('High Flow: Dead-end Channel')
