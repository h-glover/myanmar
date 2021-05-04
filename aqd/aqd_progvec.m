%% High flow tidal excursion:
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_sep17_hc.mat'),hc=aqd; hc.sal = zeros(size(hc.time));
load('aqd\aqd_sep17_lc.mat'),lc=aqd;clear aqd
load('mmi_discharge\channel_trace.mat')

% Tidal Excursion (tidal velocity amplitude * tidal period / pi
dt = (24*60*60)*(lc.time(2)-lc.time(1));
tcycle = round(12.4*60*60/dt);
numcycles=1;
s1=555;% lc = 555
e1=round(tcycle*numcycles + s1);
s2=1690;% hc = 245
e2=round(tcycle*numcycles + s2);

figure;
plot(lc.depth),hold on,plot(e1,lc.depth(e1),'ok'),plot(s1,lc.depth(s1),'ok')
title('lc hf')
figure;
plot(hc.depth),hold on,plot(e2,hc.depth(e2),'ok'),plot(s2,hc.depth(s2),'ok')
title('hc hf')

% Tidal Velocity (amplitude of depth averaged velocity)
% lc.ut = max(lc.spd_mean(s1:e1),[],'omitnan')...
%     - min(lc.spd_mean(s1:e1),[],'omitnan');
lc.ut = max(lc.spd_mean,[],'omitnan');
lc.ex=(lc.ut*(12.41*60*60))/pi/1000; % m to km

% Tidal Velocity (amplitude of depth averaged velocity)
% hc.ut = (max(hc.spd_mean(s2:e2),[],'omitnan')...
%     - min(hc.spd_mean(s2:e2),[],'omitnan'));
hc.ut = max(hc.spd_mean,[],'omitnan');
hc.ex=(hc.ut*(12.41*60*60))/pi/1000; % m to km
% compare the calcd tidal excursion to the measured tidal excursion...not
% that similar. Why??


% calculate progressive vector for dist
% define starting poisition as the aqd locations:
[hc.x,hc.y,~] = deg2utm(15.95593, 95.27682);
[lc.x,lc.y,~] = deg2utm(15.91773, 95.25484);

lc.v1_dist = lc.x + nancumsum(lc.v1_mean(s1:e1)*dt);
lc.v2_dist = lc.y + nancumsum(lc.v2_mean(s1:e1)*dt);
hc.v1_dist = hc.x + nancumsum(hc.v1_mean(s2:e2)*dt);
hc.v2_dist = hc.y + nancumsum(hc.v2_mean(s2:e2)*dt);

% calc the total distance traveled over that time: in km
v1 = [lc.v1_dist(1);lc.v1_dist];
v2 = [lc.v2_dist(1);lc.v2_dist];
lc.ex_meas = nansum(sqrt((v1(2:end)-v1(1:end-1)).^2 + (v2(2:end)-v2(1:end-1)).^2))/1000;
v1=[];v2=[];
v1 = [hc.v1_dist(1);hc.v1_dist];
v2 = [hc.v2_dist(1);hc.v2_dist];
hc.ex_meas = nansum(sqrt((v1(2:end)-v1(1:end-1)).^2 + (v2(2:end)-v2(1:end-1)).^2))/1000;

figure(10);
scatter(hc.v1_dist,hc.v2_dist,20,hc.time(s2:e2)-hc.time(s2),'d','filled')
hold on
scatter(lc.v1_dist,lc.v2_dist,20,lc.time(s1:e1)-lc.time(s1),'d','filled')
plot(ch.x,ch.y,'k.')

lc.ut_meas = lc.ex_meas/(dt*(e1-s1))*1000;
hc.ut_meas = hc.ex_meas/(dt*(e2-s2))*1000;


%% Low flow tidal excursion:
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_mar18_hc.mat'),hc=aqd; hc.sal = zeros(size(hc.time));
load('aqd\aqd_mar18_ag.mat'),ag=aqd;
load('aqd\aqd_mar18_lc.mat'),lc=aqd;clear aqd
load('mmi_discharge\channel_trace.mat')

% Tidal Excursion (tidal velocity amplitude * tidal period / pi
dt = (24*60*60)*(lc.time(2)-lc.time(1));
tcycle = round(12*60*60/dt);
numcycles=0.5;
s1=155;% lc = 155
e1=round(tcycle*numcycles + s1);
s2=300;% hc = 300
e2=round(tcycle*numcycles + s2);
s3=200;% ag = 200
e3=round(tcycle*numcycles + s3);

figure;
plot(lc.depth),hold on,plot(e1,lc.depth(e1),'ok'),plot(s1,lc.depth(s1),'ok')
title('lc lf')
figure;
plot(hc.depth),hold on,plot(e2,hc.depth(e2),'ok'),plot(s2,hc.depth(s2),'ok')
title('hc lf')
figure;
plot(ag.depth),hold on,plot(e3,ag.depth(e3),'ok'),plot(s3,ag.depth(s3),'ok')
title('ag lf')

% Tidal Velocity (amplitude of depth averaged velocity)
% lc.ut = max(lc.spd_mean(s1:e1),[],'omitnan')...
%     - min(lc.spd_mean(s1:e1),[],'omitnan');
lc.ut = max(lc.spd_mean,[],'omitnan');
lc.ex=(lc.ut*(12.41*60*60))/pi/1000; % m to km

% Tidal Velocity (amplitude of depth averaged velocity)
% hc.ut = (max(hc.spd_mean(s2:e2),[],'omitnan')...
%     - min(hc.spd_mean(s2:e2),[],'omitnan'));
hc.ut = max(hc.spd_mean,[],'omitnan');
hc.ex=(hc.ut*(12.41*60*60))/pi/1000; % m to km

% Tidal Velocity (amplitude of depth averaged velocity)
% hc.ut = (max(hc.spd_mean(s2:e2),[],'omitnan')...
%     - min(hc.spd_mean(s2:e2),[],'omitnan'));
ag.ut = max(ag.spd_mean,[],'omitnan');
ag.ex=(ag.ut*(12.41*60*60))/pi/1000; % m to km

% calculate progressive vector for dist
% define starting poisition as the aqd locations:
[hc.x,hc.y,~] = deg2utm(15.95593, 95.27682);
[lc.x,lc.y,~] = deg2utm(15.91773, 95.25484);
[ag.x,ag.y,~] = deg2utm(16.09782, 95.32051);

lc.v1_dist = lc.x + nancumsum(lc.v1_mean(s1:e1)*dt);
lc.v2_dist = lc.y + nancumsum(lc.v2_mean(s1:e1)*dt);
hc.v1_dist = hc.x + nancumsum(hc.v1_mean(s2:e2)*dt);
hc.v2_dist = hc.y + nancumsum(hc.v2_mean(s2:e2)*dt);
ag.v1_dist = ag.x + nancumsum(ag.v1_mean(s3:e3)*dt);
ag.v2_dist = ag.y + nancumsum(ag.v2_mean(s3:e3)*dt);

% calc the total distance traveled over that time: in km
v1 = [lc.v1_dist(1);lc.v1_dist];
v2 = [lc.v2_dist(1);lc.v2_dist];
lc.ex_meas = nansum(sqrt((v1(2:end)-v1(1:end-1)).^2 + (v2(2:end)-v2(1:end-1)).^2))/1000;
v1=[];v2=[];
v1 = [hc.v1_dist(1);hc.v1_dist];
v2 = [hc.v2_dist(1);hc.v2_dist];
hc.ex_meas = nansum(sqrt((v1(2:end)-v1(1:end-1)).^2 + (v2(2:end)-v2(1:end-1)).^2))/1000;
v1=[];v2=[];
v1 = [ag.v1_dist(1);ag.v1_dist];
v2 = [ag.v2_dist(1);ag.v2_dist];
ag.ex_meas = nansum(sqrt((v1(2:end)-v1(1:end-1)).^2 + (v2(2:end)-v2(1:end-1)).^2))/1000;

figure(10);
scatter(hc.v1_dist,hc.v2_dist,20,hc.time(s2:e2)-hc.time(s2),'*')
hold on
scatter(lc.v1_dist,lc.v2_dist,20,lc.time(s1:e1)-lc.time(s1),'*')
scatter(ag.v1_dist,ag.v2_dist,20,ag.time(s1:e1)-ag.time(s1),'*')
% plot(ch.x,ch.y,'k.')

lc.ut_meas = lc.ex_meas/(dt*(e1-s1))*1000;
hc.ut_meas = hc.ex_meas/(dt*(e2-s2))*1000;


%% Sept 2019 High flow tidal excursion:
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_sep19_hc.mat'),hc=aqd; hc.sal = zeros(size(hc.time));
load('aqd\aqd_sep19_ag.mat'),ag=aqd;
load('aqd\aqd_sep19_lc.mat'),lc=aqd;clear aqd
load('mmi_discharge\channel_trace.mat')

% Tidal Excursion (tidal velocity amplitude * tidal period / pi
% dt = (24*60*60)*(lc.time(2)-lc.time(1)); % 300
dt = (24*60*60)*(hc.time(2)-hc.time(1)); % 60
tcycle = round(12*60*60/dt);
numcycles=0.5;
s1=100;% lc = 100
e1=round(tcycle*numcycles + s1);
s2=2175;% hc = 2175
e2=round(tcycle*numcycles + s2);
s3=742;% ag = 742
e3=round(tcycle*numcycles + s3);

% figure;
% plot(lc.depth),hold on,plot(e1,lc.depth(e1),'ok'),plot(s1,lc.depth(s1),'ok')
% title('lc lf')
figure;
plot(hc.depth),hold on,plot(e2,hc.depth(e2),'ok'),plot(s2,hc.depth(s2),'ok')
title('hc lf')
figure;
plot(ag.depth),hold on,plot(e3,ag.depth(e3),'ok'),plot(s3,ag.depth(s3),'ok')
title('ag lf')

% % Tidal Velocity (amplitude of depth averaged velocity)
% % lc.ut = max(lc.spd_mean(s1:e1),[],'omitnan')...
% %     - min(lc.spd_mean(s1:e1),[],'omitnan');
% lc.ut = max(lc.spd_mean,[],'omitnan');
% lc.ex=(lc.ut*(12.41*60*60))/pi/1000; % m to km

% Tidal Velocity (amplitude of depth averaged velocity)
% hc.ut = (max(hc.spd_mean(s2:e2),[],'omitnan')...
%     - min(hc.spd_mean(s2:e2),[],'omitnan'));
hc.ut = max(hc.spd_mean,[],'omitnan');
hc.ex=(hc.ut*(12.41*60*60))/pi/1000; % m to km

% Tidal Velocity (amplitude of depth averaged velocity)
% hc.ut = (max(hc.spd_mean(s2:e2),[],'omitnan')...
%     - min(hc.spd_mean(s2:e2),[],'omitnan'));
ag.ut = max(ag.spd_mean,[],'omitnan');
ag.ex=(ag.ut*(12.41*60*60))/pi/1000; % m to km

% calculate progressive vector for dist
% define starting poisition as the aqd locations:
[hc.x,hc.y,~] = deg2utm(15.95593, 95.27682);
% [lc.x,lc.y,~] = deg2utm(15.91773, 95.25484);
[ag.x,ag.y,~] = deg2utm(16.09782, 95.32051);

% lc.v1_dist = lc.x + nancumsum(lc.v1_mean(s1:e1)*dt);
% lc.v2_dist = lc.y + nancumsum(lc.v2_mean(s1:e1)*dt);
hc.v1_dist = hc.x + nancumsum(hc.v1_mean(s2:e2)*dt);
hc.v2_dist = hc.y + nancumsum(hc.v2_mean(s2:e2)*dt);
ag.v1_dist = ag.x + nancumsum(ag.v1_mean(s3:e3)*dt);
ag.v2_dist = ag.y + nancumsum(ag.v2_mean(s3:e3)*dt);

% calc the total distance traveled over that time: in km
% v1 = [lc.v1_dist(1);lc.v1_dist];
% v2 = [lc.v2_dist(1);lc.v2_dist];
% lc.ex_meas = nansum(sqrt((v1(2:end)-v1(1:end-1)).^2 + (v2(2:end)-v2(1:end-1)).^2))/1000;
v1=[];v2=[];
v1 = [hc.v1_dist(1);hc.v1_dist];
v2 = [hc.v2_dist(1);hc.v2_dist];
hc.ex_meas = nansum(sqrt((v1(2:end)-v1(1:end-1)).^2 + (v2(2:end)-v2(1:end-1)).^2))/1000;
v1=[];v2=[];
v1 = [ag.v1_dist(1);ag.v1_dist];
v2 = [ag.v2_dist(1);ag.v2_dist];
ag.ex_meas = nansum(sqrt((v1(2:end)-v1(1:end-1)).^2 + (v2(2:end)-v2(1:end-1)).^2))/1000;

figure(10);
scatter(hc.v1_dist,hc.v2_dist,20,hc.time(s2:e2)-hc.time(s2),'*')
hold on
% scatter(lc.v1_dist,lc.v2_dist,20,lc.time(s1:e1)-lc.time(s1),'*')
scatter(ag.v1_dist,ag.v2_dist,20,ag.time(s1:e1)-ag.time(s1),'*')
% plot(ch.x,ch.y,'k.')

% lc.ut_meas = lc.ex_meas/(dt*(e1-s1))*1000;
% hc.ut_meas = hc.ex_meas/(dt*(e2-s2))*1000;