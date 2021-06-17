

%% wave speed in main stem from cyclone to freda
clear all,close all,clc
load('C:\GLOVER\output\myanmar\longterminst\BogaleRiverInstruments.mat')
br=BogaleRiver; clear BogaleRiver


% figure;
% plot(br.datenum,br.MeinmahlaWaterLevel,'k'),hold on
% plot(br.datenum,br.FredaWaterLevel,'r')

% read off plot = 30 minutes for 12 km
t_diff = 60*30; % in seconds
c = 11500/t_diff; %6.39 m/s

% calc wave speed from M2 phase lag: phasediff/360*12.42 (h)
% LF tidal constituents -- ut_solv(tin,uin,[],lat,cnstit)
idx = 62160:72123;
vars_mein = ut_solv(br.datenum(idx)',br.MeinmahlaDepth(idx),[],16,'auto');
vars_fred = ut_solv(br.datenum(idx)',br.FredaDepth(idx),[],16,'auto');
c_m2_lf = 11500/((345-334)/360*(12.42*60*60));
% 8.4 m/s

% solve for HF tidal constituents -- ut_solv(tin,uin,[],lat,cnstit)
idx = 132923:144000;
vars_mein = ut_solv(br.datenum(idx)',br.MeinmahlaDepth(idx),[],16,'auto');
vars_fred = ut_solv(br.datenum(idx)',br.FredaDepth(idx),[],16,'auto');
c_m2_hf = 11500/((337-327)/360*(12.42*60*60));
% 9.3 m/s

% solve for ALL tidal constituents
vars_mein = ut_solv(br.datenum',br.MeinmahlaDepth,[],16,'auto');
vars_fred = ut_solv(br.datenum',br.FredaDepth,[],16,'auto');
c_m2_all = 11500/((342-331)/360*(12.42*60*60));
% 8.4 m/s

%% wave speed in internal channel - HIGH FLOW
clear all,close all,clc
cd C:\GLOVER\output\myanmar\aqd
load('aqd_sep17_hc.mat'),hc=aqd;
load('aqd_sep17_lc.mat'),lc=aqd;

% figure;
% plot(hc.time,hc.depth,'k'),hold on,
% plot(lc.time,lc.depth,'r')
% 0.0118 = 17 mins
% 0.0500 = 72 mins
t_diff = 0.012*24*60*60; % in seconds

ext_dist = 5300;
int_dist = 14700;

c_ext = ext_dist./t_diff; % 5.1 m/s
c_int = int_dist./t_diff; % 14.2 m/s

% solve for ALL tidal constituents
% vars_hc = ut_solv(hc.time',hc.depth',[],16,'auto');
% vars_lc = ut_solv(lc.time',lc.depth',[],16,'auto');
% c_m2_all = 11500/((342-331)/360*(12.42*60*60));


%% wave speed in internal channel - LOW FLOW
clear all,close all,clc
cd C:\GLOVER\output\myanmar\aqd
load('aqd_mar18_hc.mat'),hc=aqd;
load('aqd_mar18_lc.mat'),lc=aqd;

% figure;
% plot(hc.time,hc.depth,'k'),hold on,
% plot(lc.time,lc.depth,'r')
% 0.010
t_diff = 0.01*24*60*60; % in seconds

ext_dist = 5300;
int_dist = 14700;

c_ext = ext_dist./t_diff; % 6.1 m/s
c_int = int_dist./t_diff; % 17 m/s

% solve for ALL tidal constituents
% vars_hc = ut_solv(hc.time',hc.depth',[],16,'auto');
% vars_lc = ut_solv(lc.time',lc.depth',[],16,'auto');
c_m2_all = 5700/((342-333)/360*(12.42*60*60)); % 5.1 m/s

