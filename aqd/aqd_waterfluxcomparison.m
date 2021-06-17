%% High flow - Sept 2017:
% flood is positive, ebb is negative
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_sep17_hc.mat'),hc=aqd;
load('aqd\aqd_sep17_lc.mat'),lc=aqd;clear aqd

% find sample interval
dt = (24*60*60)*(lc.time(2)-lc.time(1)); %s
tide_sec = 12.4*60*60; %length of tide in seconds

% calc instantaneous sed flux in kg/m/s, and fill missing linearly
hc.sf = hc.spd_mean.*hc.depth; % m/s * m = m3/s/m
lc.sf = lc.spd_mean.*lc.depth; % m/s * m = m3/s/m

hc.sf = fillmissing(hc.sf,'linear');
lc.sf = fillmissing(lc.sf,'linear');

% lc.sf(lc.slope<0)=(-1)*lc.sf(lc.slope<0);
% figure;plot(lc.sf)
% %%
% pick out the best single tide and find total flux during ebb and during
% flood (ignore direction - just magnitude)
s1 = 2140;
e1=round((tide_sec/dt) + s1);
tide_sf = hc.sf(s1:e1);
slp = hc.slope(s1:e1);
ebb_sum_hc = sum(tide_sf(slp<0)*dt); % m3/m over 1 tide
fld_sum_hc = sum(tide_sf(slp>0)*dt); % m3/m over 1 tide
res_sum_hc = (fld_sum_hc - ebb_sum_hc)/tide_sec; % m3/s/m

s1 = 320; %500
e1=round((tide_sec/dt) + s1);
tide_sf = lc.sf(s1:e1);
slp = lc.slope(s1:e1);
ebb_sum_lc = sum(tide_sf(slp<0)*dt); % m3/m over 1 tide
fld_sum_lc = sum(tide_sf(slp>0)*dt); % m3/m over 1 tide
res_sum_lc = (fld_sum_lc - ebb_sum_lc)/tide_sec; % m3/s/m

% total disch for comparison (pick one):
% D = 379*(10^9)/(365*24*60*60); %m3/s, mainstem annual
% D = 25880;% m3/s, mainstem in Sept from Baronas
D = 1960; %m3/s, % Bogale in Sept, spring:

% tot import:
tot_ch_width = 668; %m
tot_import = tot_ch_width*res_sum_lc; % m*m3/m/s = m3/s

ebb_transp_lc = tot_ch_width*ebb_sum_lc/tide_sec; %0.05 of total annual
fld_transp_lc = tot_ch_width*fld_sum_lc/tide_sec; %0.06 of total annual
ebb_transp_hc = 50*ebb_sum_hc/tide_sec; %0.05 of total annual
fld_transp_hc = 50*fld_sum_hc/tide_sec; %0.06 of total annual


%% Low Flow - March 2018
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_mar18_hc.mat'),hc=aqd;
load('aqd\aqd_mar18_lc.mat'),lc=aqd;clear aqd

% find sample interval
dt = (24*60*60)*(lc.time(2)-lc.time(1)); %s
tide_sec = 12.4*60*60; %length of tide in seconds

% calc instantaneous sed flux in kg/m/s, and fill missing linearly
hc.sf = hc.spd_mean.*hc.depth; % m/s * m = m3/s/m
lc.sf = lc.spd_mean.*lc.depth; % m/s * m = m3/s/m

hc.sf = fillmissing(hc.sf,'linear');
lc.sf = fillmissing(lc.sf,'linear');

% lc.sf(lc.slope<0)=(-1)*lc.sf(lc.slope<0);
% figure;plot(lc.sf)
% %%
% pick out the best single tide and find total flux during ebb and during
% flood (ignore direction - just magnitude)
s1 = 300;
e1=round((tide_sec/dt) + s1);
tide_sf = hc.sf(s1:e1);
slp = hc.slope(s1:e1);
ebb_sum_hc = sum(tide_sf(slp<0)*dt); % m3/m over 1 tide
fld_sum_hc = sum(tide_sf(slp>0)*dt); % m3/m over 1 tide
res_sum_hc = (fld_sum_hc - ebb_sum_hc)/tide_sec; % m3/s/m

s1 = 530;
e1=round((tide_sec/dt) + s1);
tide_sf = lc.sf(s1:e1);
slp = lc.slope(s1:e1);
ebb_sum_lc = sum(tide_sf(slp<0)*dt); % m3/m over 1 tide
fld_sum_lc = sum(tide_sf(slp>0)*dt); % m3/m over 1 tide
res_sum_lc = (fld_sum_lc - ebb_sum_lc)/tide_sec; % m3/s/m

% total disch for comparison (pick one):
% D = 379*(10^9)/(365*24*60*60); %m3/s, mainstem annual
% D = 2614;% m3/s, mainstem in Mar from Baronas
D = 410; %m3/s, % Bogale in Sept, spring:

% tot import:
tot_ch_width = 668; %m
tot_import = tot_ch_width*res_sum_lc; % m*m3/m/s = m3/s

ebb_transp_lc = tot_ch_width*ebb_sum_lc/tide_sec; %0.05 of total annual
fld_transp_lc = tot_ch_width*fld_sum_lc/tide_sec; %0.06 of total annual
ebb_transp_hc = 50*ebb_sum_hc/tide_sec; %0.05 of total annual
fld_transp_hc = 50*fld_sum_hc/tide_sec; %0.06 of total annual

%% High FLow - Sept 2019
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_sep19_hc.mat'),hc=aqd;
load('aqd\aqd_sep19_lc.mat'),lc=aqd;clear aqd

% find sample interval
dt_hc = (24*60*60)*(hc.time(2)-hc.time(1));
dt_lc = (24*60*60)*(lc.time(2)-lc.time(1));
tide_sec = 12.4*60*60; %length of tide in seconds

% calc instantaneous sed flux in kg/m/s, and fill missing linearly
hc.sf = hc.spd_mean.*hc.depth; % m/s * m = m3/s/m
lc.sf = lc.spd_mean.*lc.depth; % m/s * m = m3/s/m

hc.sf = fillmissing(hc.sf,'linear');
lc.sf = fillmissing(lc.sf,'linear');

% lc.sf(lc.slope<0)=(-1)*lc.sf(lc.slope<0);
% figure;plot(lc.sf)
% %%
% pick out the best single tide and find total flux during ebb and during
% flood (ignore direction - just magnitude)
s1 = 700;
e1=round((tide_sec/dt_hc) + s1);
tide_sf = hc.sf(s1:e1);
slp = hc.slope(s1:e1);
ebb_sum_hc = sum(tide_sf(slp<0)*dt_hc); % m3/m over 1 tide
fld_sum_hc = sum(tide_sf(slp>0)*dt_hc); % m3/m over 1 tide
res_sum_hc = (fld_sum_hc - ebb_sum_hc)/tide_sec; % m3/s/m

s1 = 250; 
e1=round((tide_sec/dt_lc) + s1);
tide_sf = lc.sf(s1:e1);
slp = lc.slope(s1:e1);
ebb_sum_lc = sum(tide_sf(slp<0)*dt_lc); % m3/m over 1 tide
fld_sum_lc = sum(tide_sf(slp>0)*dt_lc); % m3/m over 1 tide
res_sum_lc = (fld_sum_lc - ebb_sum_lc)/tide_sec; % m3/s/m

% total disch for comparison (pick one):
% D = 379*(10^9)/(365*24*60*60); %m3/s, mainstem annual
% D = 25880;% m3/s, mainstem in Sept from Baronas
D = 1960; %m3/s, % Bogale in Sept, spring:

% tot import:
tot_ch_width = 668; %m
tot_import = tot_ch_width*res_sum_lc; % m*m3/m/s = m3/s

ebb_transp_lc = tot_ch_width*ebb_sum_lc/tide_sec; %0.05 of total annual
fld_transp_lc = tot_ch_width*fld_sum_lc/tide_sec; %0.06 of total annual
ebb_transp_hc = 50*ebb_sum_hc/tide_sec; %0.05 of total annual
fld_transp_hc = 50*fld_sum_hc/tide_sec; %0.06 of total annual
