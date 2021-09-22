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
hc.sf = hc.spd_mean.*hc.ssc1/1000.*hc.depth; % m/s * mg/L / 1000 = (kg/m2/s)*h = kg/s/m
lc.sf = lc.spd_mean.*lc.ssc1/1000.*lc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* h = kg/s/m

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
ebb_sum_hc = sum(tide_sf(slp<0)*dt); % kg/m over 1 tide
fld_sum_hc = sum(tide_sf(slp>0)*dt); % kg/m over 1 tide
res_sum_hc = (fld_sum_hc - ebb_sum_hc)/tide_sec; % kg/s/m

s1 = 320; %500
e1=round((tide_sec/dt) + s1);
tide_sf = lc.sf(s1:e1);
slp = lc.slope(s1:e1);
ebb_sum_lc = sum(tide_sf(slp<0)*dt); % kg/m over 1 tide
fld_sum_lc = sum(tide_sf(slp>0)*dt); % kg/m over 1 tide
res_sum_lc = (fld_sum_lc - ebb_sum_lc)/tide_sec; % kg/s/m

% total disch for comparison (pick one):
% D = 326*(10^6)*1000/(365*24*60*60); %kg/s, mainstem annual
% D = 2.33*(10^6)*1000/(24*60*60);% kg/s, mainstem in Sept from Baronas
D = 0.7*1000; %kg/s, % Bogale in Sept, spring:

% tot import:
tot_ch_width = 208; % exterior channels of MMI = 668m, Int of HC channel = 208 m
tot_import = tot_ch_width*res_sum_lc; % m*kg/m/s = kg/s
% perc_retained=100*tot_import/D; %5.6% is retained in this forest
perc_ebb_transp_lc = 100*(tot_ch_width*ebb_sum_lc/tide_sec)/D; %0.05 of total annual
perc_fld_transp_lc = 100*(tot_ch_width*fld_sum_lc/tide_sec)/D; %0.06 of total annual
perc_ebb_transp_hc = 100*(50*ebb_sum_hc/tide_sec)/D; %0.05 of total annual
perc_fld_transp_hc = 100*(50*fld_sum_hc/tide_sec)/D; %0.06 of total annual

isl_area = 138250000; %m2
sed_dep_yr = tot_import*(365*24*60*60)/isl_area; %kg/s * s/yr /m2 = kg/m2/yr
accum_rate = sed_dep_yr/850*100; %kg/m2/y / 850 kg/m3 * 100 = cm/yr (BD from Cameron 2021)
% 0.02 cm/yr
% 0.006 cm/yr for lc channels off HC survey

% delta_area = 22945356200;
% 100*isl_area/delta_area; %MMI is <0.6% of delta 
%% Low Flow - March 2018
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('aqd\aqd_mar18_hc.mat'),hc=aqd;
load('aqd\aqd_mar18_lc.mat'),lc=aqd;clear aqd

% find sample interval
dt = (24*60*60)*(lc.time(2)-lc.time(1)); %s
tide_sec = 12.4*60*60; %length of tide in seconds

% calc instantaneous sed flux in kg/m/s, and fill missing linearly
hc.sf = hc.spd_mean.*hc.ssc1/1000.*hc.depth; % m/s * mg/L / 1000 = (kg/m2/s)*h = kg/s/m
lc.sf = lc.spd_mean.*lc.ssc1/1000.*lc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* h = kg/s/m

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
ebb_sum_hc = sum(tide_sf(slp<0)*dt); % kg/m over 1 tide
fld_sum_hc = sum(tide_sf(slp>0)*dt); % kg/m over 1 tide
res_sum_hc = (fld_sum_hc - ebb_sum_hc)/tide_sec; % kg/s/m

s1 = 530; %500
e1=round((tide_sec/dt) + s1);
tide_sf = lc.sf(s1:e1);
slp = lc.slope(s1:e1);
ebb_sum_lc = sum(tide_sf(slp<0)*dt); % kg/m over 1 tide
fld_sum_lc = sum(tide_sf(slp>0)*dt); % kg/m over 1 tide
res_sum_lc = (fld_sum_lc - ebb_sum_lc)/tide_sec; % kg/s/m

% total disch for comparison (pick one):
% D = 326*(10^6)*1000/(365*24*60*60); %kg/s, mainstem annual
% D = 2.33*(10^6)*1000/(24*60*60);% kg/s, mainstem in Sept from Baronas
D = 0.1*1000; %kg/s, % Bogale in Sept, spring:

% tot import:
tot_ch_width = 208;% exterior channels of MMI = 668m, Int of HC channel = 208 m
tot_import = tot_ch_width*res_sum_lc; % m*kg/m/s = kg/s
% perc_retained=100*tot_import/D; %5.6% is retained in this forest
perc_ebb_transp_lc = 100*(tot_ch_width*ebb_sum_lc/tide_sec)/D; %0.05 of total annual
perc_fld_transp_lc = 100*(tot_ch_width*fld_sum_lc/tide_sec)/D; %0.06 of total annual
perc_ebb_transp_hc = 100*(50*ebb_sum_hc/tide_sec)/D; %0.05 of total annual
perc_fld_transp_hc = 100*(50*fld_sum_hc/tide_sec)/D; %0.06 of total annual

isl_area = 138250000; %m2
sed_dep_yr = tot_import*(365*24*60*60)/isl_area; %kg/s * s/yr /m2 = kg/m2/yr
accum_rate = sed_dep_yr/850*100; %kg/m2/y / 850 kg/m3 * 100 = cm/yr (BD from Cameron 2021)
% 0.06 cm/yr for ext LCs
% 0.02 cm/yr for Int LCs off of HC

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
hc.sf = hc.spd_mean.*hc.ssc1/1000.*hc.depth; % m/s * mg/L / 1000 = (kg/m2/s)*h = kg/s/m
lc.sf = lc.spd_mean.*lc.ssc1/1000.*lc.depth; % m/s * mg/L / 1000 = (kg/m2/s)* h = kg/s/m

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
ebb_sum_hc = sum(tide_sf(slp<0)*dt_hc); % kg/m over 1 tide
fld_sum_hc = sum(tide_sf(slp>0)*dt_hc); % kg/m over 1 tide
res_sum_hc = (fld_sum_hc - ebb_sum_hc)/tide_sec; % kg/s/m

s1 = 250;
e1=round((tide_sec/dt_lc) + s1);
tide_sf = lc.sf(s1:e1);
slp = lc.slope(s1:e1);
ebb_sum_lc = sum(tide_sf(slp<0)*dt_lc); % kg/m over 1 tide
fld_sum_lc = sum(tide_sf(slp>0)*dt_lc); % kg/m over 1 tide
res_sum_lc = (fld_sum_lc - ebb_sum_lc)/tide_sec; % kg/s/m

% total disch for comparison (pick one):
% D = 326*(10^6)*1000/(365*24*60*60); %kg/s, mainstem annual
% D = 2.33*(10^6)*1000/(24*60*60);% kg/s, mainstem in Sept from Baronas
D = 0.7*1000; %kg/s, % Bogale in Sept, spring:

% tot import:
tot_ch_width = 208;% exterior channels of MMI = 668m, Int of HC channel = 208 m
tot_import = tot_ch_width*res_sum_lc; % m*kg/m/s = kg/s
% perc_retained=100*tot_import/D; %5.6% is retained in this forest
perc_ebb_transp_lc = 100*(tot_ch_width*ebb_sum_lc/tide_sec)/D; %0.05 of total annual
perc_fld_transp_lc = 100*(tot_ch_width*fld_sum_lc/tide_sec)/D; %0.06 of total annual
perc_ebb_transp_hc = 100*(50*ebb_sum_hc/tide_sec)/D; %0.05 of total annual
perc_fld_transp_hc = 100*(50*fld_sum_hc/tide_sec)/D; %0.06 of total annual

isl_area = 138250000; %m2
sed_dep_yr = tot_import*(365*24*60*60)/isl_area; %kg/s * s/yr /m2 = kg/m2/yr
accum_rate = sed_dep_yr/850*100; %kg/m2/y / 850 kg/m3 * 100 = cm/yr (BD from Cameron 2021)

% 0.12 cm/yr
% 0 cm/yr for Int LCs off of HC