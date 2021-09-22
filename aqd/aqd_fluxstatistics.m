%% SED FLUX (WATER STATS ARE BELOW)
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

hc.t_idx = 1:floor(tide_sec/dt):length(hc.sf);
lc.t_idx = 1:floor(tide_sec/dt):length(lc.sf);

for jj=1:length(hc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = hc.t_idx(jj);e1 = hc.t_idx(jj+1)-1;
    tide_sf = hc.sf(s1:e1);
    slp = hc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    hc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt)/tide_sec; % kg/m over 1 tide
    hc.fld_sum(jj) = sum(tide_sf(slp>0)*dt)/tide_sec; % kg/m over 1 tide
end

hc.res_sum = (hc.fld_sum - hc.ebb_sum); % kg/m (divide by seconds in tide to get kg/m/s

t(1)=nanmean(hc.ebb_sum);
t(2)=std(hc.ebb_sum,'omitnan');
t(3)=nanmean(hc.fld_sum);
t(4)=std(hc.fld_sum,'omitnan');
t(5)=nanmean(hc.res_sum);
t(6)=std(hc.res_sum,'omitnan');
nanmean(hc.res_sum)*50/700*100

for jj=1:length(lc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = lc.t_idx(jj);e1 = lc.t_idx(jj+1)-1;
    tide_sf = lc.sf(s1:e1);
    slp = lc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    lc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt)/tide_sec; % kg/m over 1 tide
    lc.fld_sum(jj) = sum(tide_sf(slp>0)*dt)/tide_sec; % kg/m over 1 tide
end
lc.res_sum = (lc.fld_sum - lc.ebb_sum); % kg/m/s

tl(1)=nanmean(lc.ebb_sum);
tl(2)=std(lc.ebb_sum,'omitnan');
tl(3)=nanmean(lc.fld_sum);
tl(4)=std(lc.fld_sum,'omitnan');
tl(5)=nanmean(lc.res_sum);
tl(6)=std(lc.res_sum,'omitnan');

nanmean(lc.res_sum)*668
% nanmean(lc.res_sum)*668 
% nanmean(lc.res_sum)*208 
%% low flow: march 2018
% flood is positive, ebb is negative
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

hc.t_idx = 1:floor(tide_sec/dt):length(hc.sf);
lc.t_idx = 1:floor(tide_sec/dt):length(lc.sf);

for jj=1:length(hc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = hc.t_idx(jj);e1 = hc.t_idx(jj+1)-1;
    tide_sf = hc.sf(s1:e1);
    slp = hc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    hc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt)/tide_sec; % kg/m over 1 tide
    hc.fld_sum(jj) = sum(tide_sf(slp>0)*dt)/tide_sec; % kg/m over 1 tide
end

hc.res_sum = (hc.fld_sum - hc.ebb_sum); % kg/m (divide by seconds in tide to get kg/m/s

t(1)=nanmean(hc.ebb_sum);
t(2)=std(hc.ebb_sum,'omitnan');
t(3)=nanmean(hc.fld_sum);
t(4)=std(hc.fld_sum,'omitnan');
t(5)=nanmean(hc.res_sum);
t(6)=std(hc.res_sum,'omitnan');
nanmean(hc.res_sum)*50/120*100

for jj=1:length(lc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = lc.t_idx(jj);e1 = lc.t_idx(jj+1)-1;
    tide_sf = lc.sf(s1:e1);
    slp = lc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    lc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt)/tide_sec; % kg/m over 1 tide
    lc.fld_sum(jj) = sum(tide_sf(slp>0)*dt)/tide_sec; % kg/m over 1 tide
end
lc.res_sum = (lc.fld_sum - lc.ebb_sum); % kg/m/s

tl(1)=nanmean(lc.ebb_sum);
tl(2)=std(lc.ebb_sum,'omitnan');
tl(3)=nanmean(lc.fld_sum);
tl(4)=std(lc.fld_sum,'omitnan');
tl(5)=nanmean(lc.res_sum);
tl(6)=std(lc.res_sum,'omitnan');

nanmean(lc.res_sum)*668%/120*100

%% sept 2019
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

hc.t_idx = 1:floor(tide_sec/dt_hc):length(hc.sf);
lc.t_idx = 1:floor(tide_sec/dt_lc):length(lc.sf);

for jj=1:length(hc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = hc.t_idx(jj);e1 = hc.t_idx(jj+1)-1;
    tide_sf = hc.sf(s1:e1);
    slp = hc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    hc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt_hc)/tide_sec; % kg/m over 1 tide
    hc.fld_sum(jj) = sum(tide_sf(slp>0)*dt_hc)/tide_sec; % kg/m over 1 tide
end

hc.res_sum = (hc.fld_sum - hc.ebb_sum); % kg/m (divide by seconds in tide to get kg/m/s

t(1)=nanmean(hc.ebb_sum);
t(2)=std(hc.ebb_sum,'omitnan');
t(3)=nanmean(hc.fld_sum);
t(4)=std(hc.fld_sum,'omitnan');
t(5)=nanmean(hc.res_sum);
t(6)=std(hc.res_sum,'omitnan');
nanmean(hc.res_sum)*50/700*100

for jj=1:length(lc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = lc.t_idx(jj);e1 = lc.t_idx(jj+1)-1;
    tide_sf = lc.sf(s1:e1);
    slp = lc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    lc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt_lc)/tide_sec; % kg/m over 1 tide
    lc.fld_sum(jj) = sum(tide_sf(slp>0)*dt_lc)/tide_sec; % kg/m over 1 tide
end
lc.res_sum = (lc.fld_sum - lc.ebb_sum); % kg/m/s

tl(1)=nanmean(lc.ebb_sum);
tl(2)=std(lc.ebb_sum,'omitnan');
tl(3)=nanmean(lc.fld_sum);
tl(4)=std(lc.fld_sum,'omitnan');
tl(5)=nanmean(lc.res_sum);
tl(6)=std(lc.res_sum,'omitnan');

nanmean(lc.res_sum)*668%/700*100

%%
% WATER FLUX

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

hc.t_idx = 1:floor(tide_sec/dt):length(hc.sf);
lc.t_idx = 1:floor(tide_sec/dt):length(lc.sf);

for jj=1:length(hc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = hc.t_idx(jj);e1 = hc.t_idx(jj+1)-1;
    tide_sf = hc.sf(s1:e1);
    slp = hc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    hc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt)/tide_sec; % kg/m over 1 tide
    hc.fld_sum(jj) = sum(tide_sf(slp>0)*dt)/tide_sec; % kg/m over 1 tide
end

t(1)=nanmean(hc.ebb_sum);
t(2)=std(hc.ebb_sum,'omitnan');
t(3)=nanmean(hc.fld_sum);
t(4)=std(hc.fld_sum,'omitnan');

t(5)=(t(1)*50)/1960*100;
t(6)=(t(3)*50)/1960*100;



for jj=1:length(lc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = lc.t_idx(jj);e1 = lc.t_idx(jj+1)-1;
    tide_sf = lc.sf(s1:e1);
    slp = lc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    lc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt)/tide_sec; % kg/m over 1 tide
    lc.fld_sum(jj) = sum(tide_sf(slp>0)*dt)/tide_sec; % kg/m over 1 tide
end


tl(1)=nanmean(lc.ebb_sum);
tl(2)=std(lc.ebb_sum,'omitnan');
tl(3)=nanmean(lc.fld_sum);
tl(4)=std(lc.fld_sum,'omitnan');

tl(5)=(tl(1)*668)/1960*100;
tl(6)=(tl(3)*668)/1960*100;


%% Low flow - Mar 2018:
% flood is positive, ebb is negative
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

hc.t_idx = 1:floor(tide_sec/dt):length(hc.sf);
lc.t_idx = 1:floor(tide_sec/dt):length(lc.sf);

for jj=1:length(hc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = hc.t_idx(jj);e1 = hc.t_idx(jj+1)-1;
    tide_sf = hc.sf(s1:e1);
    slp = hc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    hc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt)/tide_sec; % kg/m over 1 tide
    hc.fld_sum(jj) = sum(tide_sf(slp>0)*dt)/tide_sec; % kg/m over 1 tide
end

t(1)=nanmean(hc.ebb_sum);
t(2)=std(hc.ebb_sum,'omitnan');
t(3)=nanmean(hc.fld_sum);
t(4)=std(hc.fld_sum,'omitnan');

t(5)=(t(1)*50)/410*100;
t(6)=(t(3)*50)/410*100;



for jj=1:length(lc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = lc.t_idx(jj);e1 = lc.t_idx(jj+1)-1;
    tide_sf = lc.sf(s1:e1);
    slp = lc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    lc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt)/tide_sec; % kg/m over 1 tide
    lc.fld_sum(jj) = sum(tide_sf(slp>0)*dt)/tide_sec; % kg/m over 1 tide
end


tl(1)=nanmean(lc.ebb_sum);
tl(2)=std(lc.ebb_sum,'omitnan');
tl(3)=nanmean(lc.fld_sum);
tl(4)=std(lc.fld_sum,'omitnan');

tl(5)=(tl(1)*668)/410*100;
tl(6)=(tl(3)*668)/410*100;

%% Low flow - Sept 2019:
% flood is positive, ebb is negative
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

hc.t_idx = 1:floor(tide_sec/dt_hc):length(hc.sf);
lc.t_idx = 1:floor(tide_sec/dt_lc):length(lc.sf);

for jj=1:length(hc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = hc.t_idx(jj);e1 = hc.t_idx(jj+1)-1;
    tide_sf = hc.sf(s1:e1);
    slp = hc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    hc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt_hc)/tide_sec; % kg/m over 1 tide
    hc.fld_sum(jj) = sum(tide_sf(slp>0)*dt_hc)/tide_sec; % kg/m over 1 tide
end

t(1)=nanmean(hc.ebb_sum);
t(2)=std(hc.ebb_sum,'omitnan');
t(3)=nanmean(hc.fld_sum);
t(4)=std(hc.fld_sum,'omitnan');

t(5)=(t(1)*50)/1960*100;
t(6)=(t(3)*50)/1960*100;



for jj=1:length(lc.t_idx)-1
    % pull out the chuck corresponding to that ebb/flood cycle
    s1 = lc.t_idx(jj);e1 = lc.t_idx(jj+1)-1;
    tide_sf = lc.sf(s1:e1);
    slp = lc.slope(s1:e1);
    
    % calc total flux per ebb and per flood in kg/m
    lc.ebb_sum(jj) = sum(tide_sf(slp<0)*dt_lc)/tide_sec; % kg/m over 1 tide
    lc.fld_sum(jj) = sum(tide_sf(slp>0)*dt_lc)/tide_sec; % kg/m over 1 tide
end


tl(1)=nanmean(lc.ebb_sum);
tl(2)=std(lc.ebb_sum,'omitnan');
tl(3)=nanmean(lc.fld_sum);
tl(4)=std(lc.fld_sum,'omitnan');

tl(5)=(tl(1)*668)/1960*100;
tl(6)=(tl(3)*668)/1960*100;



