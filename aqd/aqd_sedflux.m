%% High flow tidal excursion:
% flood is positive, ebb is negative
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_sep17_hc.mat'),hc=aqd; hc.sal = zeros(size(hc.time));
load('aqd\aqd_sep17_lc.mat'),lc=aqd;clear aqd

dt = (24*60*60)*(lc.time(2)-lc.time(1));

hc_vel = hc.spd_mean;
hc_vel(hc.slope<0)=(-1)*hc_vel(hc.slope<0);

lc_vel = lc.spd_mean;
lc_vel(lc.slope<0)=(-1)*lc_vel(lc.slope<0);

hc.sf = hc_vel.*hc.ssc1/1000.*hc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* (1m*dh) = kg/s
lc.sf = lc_vel.*lc.ssc1/1000.*lc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* (1m*dh) = kg/s

hc.sf = fillmissing(hc.sf,'linear');
lc.sf = fillmissing(lc.sf,'linear');

hc.sf_cum = nancumsum(hc.sf*dt);
lc.sf_cum = nancumsum(lc.sf*dt);
numtides = (lc.time(end)-lc.time(1))/(12.4/24);


figure;
subplot(211)
yyaxis left,plot(hc.depth),hold on
yyaxis right,plot(hc.sf_cum),hold on
subplot(212)
yyaxis left,plot(lc.depth)
yyaxis right,plot(lc.sf_cum)
figure;
subplot(211)
yyaxis left,plot(hc.depth),hold on
yyaxis right,plot(hc.sf),hold on
subplot(212)
yyaxis left,plot(lc.depth)
yyaxis right,plot(lc.sf)
% 
% s1 = 2200;
% e1=round((1*12.4*60*60/dt) + s1);
% test = nancumsum(hc.sf(s1:e1)*dt);
% figure;
% subplot(411),plot(hc.depth(s1:e1))
% subplot(412),plot(hc_vel(s1:e1))
% subplot(413),plot(hc.ssc1(s1:e1))
% subplot(414),plot(test)
% 
% figure;
% subplot(411),plot(lc.depth)
% subplot(412),plot(lc_vel)
% subplot(413),plot(lc.ssc1)
% subplot(414),plot(lc.sf)
% 
% s1 = 500;
% e1=round((5*12.4*60*60/dt) + s1);
% test = nancumsum(lc.sf(s1:e1)*dt);
% figure;
% subplot(411),plot(lc.depth(s1:e1))
% subplot(412),plot(lc_vel(s1:e1))
% subplot(413),plot(lc.ssc1(s1:e1))
% subplot(414),plot(test)

% calc the residual tidal flux as the int avg over 1 tidal cycle:
idx = 2000:2000+((12.4*60*60)/dt);
hc.res = nanmean(hc.sf(idx));
lc.res = nanmean(lc.sf(idx));

%% Low flow tidal excursion:
% flood is positive, ebb is negative
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_mar18_hc.mat'),hc=aqd; hc.sal = zeros(size(hc.time));
load('aqd\aqd_mar18_ag.mat'),ag=aqd;
load('aqd\aqd_mar18_lc.mat'),lc=aqd;clear aqd

dt = (24*60*60)*(lc.time(2)-lc.time(1));

hc_vel = hc.spd_mean;
hc_vel(hc.slope<0)=(-1)*hc_vel(hc.slope<0);

lc_vel = lc.spd_mean;
lc_vel(lc.slope<0)=(-1)*lc_vel(lc.slope<0);

ag_vel = ag.spd_mean;
ag_vel(ag.slope<0)=(-1)*ag_vel(ag.slope<0);

hc.sf = hc_vel.*hc.ssc1/1000.*hc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* (1m*dh) = kg/s
lc.sf = lc_vel.*lc.ssc1/1000.*lc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* (1m*dh) = kg/s
ag.sf = ag_vel.*ag.ssc1/1000.*ag.depth;

hc.sf = fillmissing(hc.sf,'linear');
lc.sf = fillmissing(lc.sf,'linear');
ag.sf = fillmissing(ag.sf,'linear');

hc.sf_cum = nancumsum(hc.sf*dt);
lc.sf_cum = nancumsum(lc.sf*dt);
ag.sf_cum = nancumsum(ag.sf*dt);

figure;
subplot(311)
yyaxis left,plot(hc.depth),title('hc')
yyaxis right,plot(hc.sf_cum)
subplot(312)
yyaxis left,plot(lc.depth),title('lc')
yyaxis right,plot(lc.sf_cum)
subplot(313)
yyaxis left,plot(ag.depth),title('ag')
yyaxis right,plot(ag.sf_cum)
figure;
subplot(311)
yyaxis left,plot(hc.depth),title('hc')
yyaxis right,plot(hc.sf)
subplot(312)
yyaxis left,plot(lc.depth),title('lc')
yyaxis right,plot(lc.sf)
subplot(313)
yyaxis left,plot(ag.depth),title('ag')
yyaxis right,plot(ag.sf)

figure;
subplot(411),plot(lc.depth),title('lc')
subplot(412),plot(lc_vel)
subplot(413),plot(lc.ssc1,'k'),hold on,plot(lc.ssc2,'r')
subplot(414),plot(lc.sf)
figure;
subplot(411),plot(hc.depth),title('hc')
subplot(412),plot(hc_vel)
subplot(413),plot(hc.ssc1,'k'),hold on,plot(hc.ssc2,'r')
subplot(414),plot(hc.sf)
figure;
subplot(411),plot(ag.depth),title('ag')
subplot(412),plot(ag_vel)
subplot(413),plot(ag.ssc1,'k'),hold on,plot(ag.ssc2,'r')
subplot(414),plot(ag.sf)


idx = round(3*12.4*60*60/dt);
test_avg=ag.sf_cum(idx)/3;


% % calc the residual tidal flux as the int avg over 1 tidal cycle:
idx = 300:300+((12.4*60*60)/dt);
hc.res = nanmean(hc.sf(idx));
idx = 600:600+((12.4*60*60)/dt);
lc.res = nanmean(lc.sf(idx));
idx = 100:100+((12.4*60*60)/dt);
ag.res = nanmean(ag.sf(idx));



%% second high flow tidal excursion:
% flood is positive, ebb is negative
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_sep19_hc.mat'),hc=aqd; hc.sal = zeros(size(hc.time));
load('aqd\aqd_sep19_ag.mat'),ag=aqd;
load('aqd\aqd_sep19_lc.mat'),lc=aqd;clear aqd

dt = (24*60*60)*(hc.time(2)-hc.time(1));
dt_lc = (24*60*60)*(lc.time(2)-lc.time(1));

hc_vel = hc.spd_mean;
hc_vel(hc.slope<0)=(-1)*hc_vel(hc.slope<0);

lc_vel = lc.spd_mean;
lc_vel(lc.slope<0)=(-1)*lc_vel(lc.slope<0);

ag_vel = ag.spd_mean;%ag_vel(ag_vel>0.7)=NaN;
ag_vel(ag.slope<0)=(-1)*ag_vel(ag.slope<0);

hc.sf = hc_vel.*hc.ssc2/1000.*hc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* (1m*dh) = kg/s
lc.sf = lc_vel.*lc.ssc2/1000.*lc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* (1m*dh) = kg/s
ag.sf = ag_vel.*ag.ssc1/1000.*ag.depth;

hc.sf = fillmissing(hc.sf,'linear');
lc.sf = fillmissing(lc.sf,'linear');
% ag.sf = fillmissing(ag.sf,'linear');

hc.sf_cum = nancumsum(hc.sf*dt);
lc.sf_cum = nancumsum(lc.sf*dt_lc);
ag.sf_cum = nancumsum(ag.sf*dt);

figure;
subplot(311)
yyaxis left,plot(hc.depth),title('hc')
yyaxis right,plot(hc.sf_cum)
subplot(312)
yyaxis left,plot(lc.depth),title('lc')
yyaxis right,plot(lc.sf_cum)
subplot(313)
yyaxis left,plot(ag.depth),title('ag')
yyaxis right,plot(ag.sf_cum)
figure;
subplot(311)
yyaxis left,plot(hc.depth),title('hc')
yyaxis right,plot(hc.sf)
subplot(312)
yyaxis left,plot(lc.depth),title('lc')
yyaxis right,plot(lc.sf)
subplot(313)
yyaxis left,plot(ag.depth),title('ag')
yyaxis right,plot(ag.sf)

figure;
subplot(411),plot(lc.depth),title('lc')
subplot(412),plot(lc_vel)
subplot(413),plot(lc.ssc1,'k'),hold on,plot(lc.ssc2,'r')
subplot(414),plot(lc.sf)
figure;
subplot(411),plot(hc.depth),title('hc')
subplot(412),plot(hc_vel)
subplot(413),plot(hc.ssc1,'k'),hold on,plot(hc.ssc2,'r')
subplot(414),plot(hc.sf)
figure;
subplot(411),plot(ag.depth),title('ag')
subplot(412),plot(ag_vel)
subplot(413),plot(ag.ssc1,'k'),hold on,plot(ag.ssc2,'r')
subplot(414),plot(ag.sf)


% idx = round(3*12.4*60*60/dt_lc);
% lc.sf_cum(idx)/3


% calc the residual tidal flux as the int avg over 1 tidal cycle:
idx = 700:700+((12.4*60*60)/dt);
hc.res = nanmean(hc.sf(idx));
idx = 250:250+((12.4*60*60)/dt_lc);
lc.res = nanmean(lc.sf(idx));
idx = 1:1+((12.4*60*60)/dt);
ag.res = nanmean(ag.sf(idx));