% plot of data overlaps
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('mmi_discharge\br_int_march2018.mat')
load('aqd\aqd_mar18_hc.mat'),hc=aqd;
load('aqd\aqd_mar18_lc.mat'),lc=aqd;clear aqd
%%
hc.bindepth = hc.depth-0.5 - repmat(hc.z,[length(hc.time),1]);
% hc.bindepth(hc.bindepth<0)=NaN;
lc.bindepth = lc.depth-0.5 - repmat(lc.z,[length(lc.time),1]);
% lc.bindepth(lc.bindepth<0)=NaN;
hc.t = repmat(hc.time,[1,20]);
lc.t = repmat(lc.time,[1,20]);

figure;
plot(hc.time,hc.depth,'k'),hold on
plot(lc.time,lc.depth,'r')
yyaxis right
for jj=1:3
    plot(adcp(jj).time,adcp(jj).elapdist)
end
legend({'hc','lc','a1','a2','a3'})
% adcp(3) matches with the HC channel aqd
% adcp(1,2) matches with the LC channel aqd
%%
figure;
subplot(221),pcolor(adcp(3).time,adcp(3).z(:,1),adcp(3).spd./100)
ylim([0 10]),caxis([0 0.5])
subplot(222),pcolor(adcp(3).time,adcp(3).z(:,1),adcp(3).dir)
ylim([0 10]),caxis([0 360])
subplot(223),pcolor(hc.t,hc.bindepth,hc.spd),caxis([0 0.5])
subplot(224),pcolor(hc.t,hc.bindepth,hc.dir),caxis([0 360])

for jj=1:4
    subplot(2,2,jj),shading flat,colorbar,axis ij
    xlim([adcp(3).time(1) adcp(3).time(end)])
    datetick('x','HH:MM','keeplimits')
end

%% TS diagram:

figure;
scatter(hc.sal,hc.temp,'*'),hold on
scatter(lc.sal,lc.temp,'o')


figure;
scatter(adcp(1).sal,adcp(1).theta,'*'),hold on
scatter(adcp(3).sal,adcp(3).theta,'*')