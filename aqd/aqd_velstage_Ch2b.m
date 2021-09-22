clear all,close all,clc
cd C:\GLOVER\output\myanmar\aqd

load('aqd_sep17_lc.mat')
aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.05;
aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
lc17 = aqd;clear aqd

load('aqd_mar18_lc.mat')
aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.08;
aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
lc18=aqd;clear aqd

load('aqd_sep19_ag.mat')
aqd.spd_mean(aqd.spd_mean>0.3)=NaN;
aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.07;
aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
aqd.ssc1=aqd.ssc1*0.2;
ag19 = aqd;clear aqd

load('aqd_mar18_ag.mat')
aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.05;
aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
ag18=aqd;clear aqd



figure; 
subplot(221)
scatter(lc17.spd_mean,lc17.depth,[],lc17.ssc1,'.')
title('HF MMI')
subplot(223)
scatter(lc18.spd_mean,lc18.depth,[],lc18.ssc1,'.')
title('LF MMI')

subplot(222)
idx = 1:length(ag19.time);
scatter(ag19.spd_mean(idx),ag19.depth(idx),[],ag19.ssc1(idx),'.')
title('HF AG')
subplot(224)
idx = 1:length(ag18.time)-80;
scatter(ag18.spd_mean(idx),ag18.depth(idx),[],ag18.ssc1(idx),'.')
title('LF AG')
colorbar

for jj=1:4
    subplot(2,2,jj)
    caxis([0 100]),colormap(cmocean('turbid'))
    axis([-0.5 0.5 0 3.5])
    ylabel('Water Depth (m)')
    xlabel('Water Velocity (m/s)')
end
