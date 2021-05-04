clear all,close all,clc

cd C:\GLOVER\output\myanmar\longterminst

load('BogaleRiverWell'),bw=BogaleWell;
load('BogaleRiverInstruments'),br=BogaleRiver;
load('BogaleRiverWeather')

bw.MeinmahlaWellTemp=movmean(bw.MeinmahlaWellTemp,7);


figure;
subplot(2,2,1:2)
plot(br.datenum,br.MeinmahlaWaterLevel,'b'),hold on
plot(bw.datenum,bw.MeinmahlaWellDepth-0.35,'k--')
ylim([-2 2]),ylabel('water level')
yyaxis right
plot(bw.datenum,bw.MeinmahlaWellTemp,'r')
ylim([25 29]),ylabel('temp (C)')
xlim([datenum('3/1/2017') datenum('10/1/2017')])
datetick('x','mm yy','keeplimits')
legend({'River','Well','Well Temp'})

subplot(2,2,3)
plot(br.datenum,br.MeinmahlaWaterLevel,'b'),hold on
plot(bw.datenum,bw.MeinmahlaWellDepth-0.35,'k')
ylim([-2 1.5])
yyaxis right
plot(bw.datenum,bw.MeinmahlaWellTemp,'r')
ylim([25 29])
xlim([bw.datenum(7875) bw.datenum(10700)])
datetick('x','mm dd','keeplimits')
title('dry season')

subplot(2,2,4)
plot(br.datenum,br.MeinmahlaWaterLevel,'b'),hold on
plot(bw.datenum,bw.MeinmahlaWellDepth-0.35,'k')
ylim([-2 1.5])
yyaxis right
plot(bw.datenum,bw.MeinmahlaWellTemp,'r')
ylim([25 29])
xlim([bw.datenum(22100) bw.datenum(24925)])
datetick('x','mm dd','keeplimits')
title('wet season')

%% fig for s2s talk:
clear all,close all,clc

cd C:\GLOVER\output\myanmar\longterminst

load('BogaleRiverWell'),bw=BogaleWell;
load('BogaleRiverInstruments'),br=BogaleRiver;
load('BogaleRiverWeather')

bw.MeinmahlaWellTemp=movmean(bw.MeinmahlaWellTemp,7);

figure;
subplot(1,2,1)
plot(br.datenum,br.MeinmahlaWaterLevel,'b'),hold on
plot(bw.datenum,bw.MeinmahlaWellDepth-0.35,'k')
ylim([-2 1.5])
ylabel('Water level (m)')
xlim([bw.datenum(7875) bw.datenum(10700)])
datetick('x','mmm dd','keeplimits')
title('LF (NE) Monsoon')

subplot(1,2,2)
plot(br.datenum,br.MeinmahlaWaterLevel,'b'),hold on
plot(bw.datenum,bw.MeinmahlaWellDepth-0.35,'k')
ylim([-2 1.5])
xlim([bw.datenum(22100) bw.datenum(24925)])
datetick('x','mmm dd','keeplimits')
ylabel('Water level (m)')
title('HF (SW) Monsoon')
legend({'River','Water table'})