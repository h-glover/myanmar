%% calc sed flux at CTD sites: 
% find vel within 30m of ctd profile
% avg vel profile
% multiply by ssc
% how to get cross-sectional area of channel? or per/m width?

clear all,close all,clc

cd C:\GLOVER\output\myanmar
load('mmi_discharge\br_int_march2018.mat')
load('ctdcasts\AyeMar18_CTD_all.mat')

for jj=1:length(CombinedProfiles)
    alltime_ctd(jj)=CombinedProfiles(jj).time(1);
end
[alltime_ctd,idx]=sort(alltime_ctd,'ascend');
ctd=CombinedProfiles(idx);

% concatenate the profile indexs and times
allidx_ctd=1:length(ctd);
alltime_ctd=datevec(alltime_ctd);alltime_ctd(:,6)=0;alltime_ctd=datenum(alltime_ctd);

% assign a transect number to each time point
for jj=1:length(adcp)
    adcp(jj).transidx=jj*ones(1,length(adcp(jj).time));
    adcp(jj).depth=nanmean(adcp(jj).depth);
    adcp(jj).depth(adcp(jj).depth<0)=NaN;
    adcp(jj).depth=fillmissing(adcp(jj).depth,'nearest');
end

% concatenate into 1 vector and remove seconds
allidx_adcp=horzcat(adcp.transidx);
alltime_adcp=horzcat(adcp.time);
alltime_adcp=datevec(alltime_adcp);alltime_adcp(:,6)=0;
alltime_adcp=datenum(alltime_adcp);

% compare adcp time and ctd time (seconds removed)
[~,iA,iC]=intersect(alltime_adcp,alltime_ctd);

% pull out all ctd profiles from this day
ctd_trans=allidx_ctd(iC);
allidx_adcp=allidx_adcp(iA);
ctd=ctd(iC);
alltime_ctd=alltime_ctd(iC);

% interp to a uniform depth interval
intrp_depth =0:0.05:15; %max adcp depth is 15m

for jj=1:length(ctd)
    % interp ssc and salinity into the new uniform depth matrix
    ctd(jj).SSCCal = interp1(ctd(jj).Depth,ctd(jj).SSCCal,intrp_depth)';
    ctd(jj).Salinity = interp1(ctd(jj).Depth,ctd(jj).Salinity,intrp_depth)';
    ctd(jj).SSCCal(ctd(jj).SSCCal>1500)=1000;
    ctd(jj).Temperature = interp1(ctd(jj).Depth,ctd(jj).Temperature,intrp_depth)';

    % get the adcp transect
    ctd(jj).transect=allidx_adcp(jj);
    L(jj)=length(ctd(jj).time);
end
% % brute force...
% ctd1=ctd([1:7]);
% ctd3=ctd([8:end]);

% first map ctd1 onto adcp transect1
for kk=[1,3]
adcp(kk).ssc = NaN(length(intrp_depth),length(adcp(kk).time));
adcp(kk).sal = NaN(length(intrp_depth),length(adcp(kk).time));
adcp(kk).theta = NaN(length(intrp_depth),length(adcp(kk).time));
adcp(kk).sed_flux = NaN(size(adcp(kk).east));
adcp(kk).sal_flux = NaN(size(adcp(kk).east));

[~,iC,iA] = intersect(alltime_ctd,adcp(kk).time);
for jj=1:length(iC)
    adcp(kk).ssc(:,iA(jj)) = ctd(jj).SSCCal;
    adcp(kk).sal(:,iA(jj)) = ctd(jj).Salinity;
    adcp(kk).theta(:,iA(jj)) = ctd(jj).Temperature;
    % calculate sed flux - first put in ssc, later multiply by vel components
    adcp(kk).sed_flux(:,iA(jj)) = interp1(intrp_depth,adcp(kk).ssc(:,iA(jj)),adcp(kk).z(:,1))./1000; %convert to kg/m3
    adcp(kk).sal_flux(:,iA(jj)) = interp1(intrp_depth,adcp(kk).sal(:,iA(jj)),adcp(kk).z(:,1));
end

% interpolate to fill the grids
adcp(kk).sal = fillmissing(adcp(kk).sal,'nearest',1);
adcp(kk).sal = fillmissing(adcp(kk).sal,'linear',2);
adcp(kk).ssc = fillmissing(adcp(kk).ssc,'nearest',1);
adcp(kk).ssc = fillmissing(adcp(kk).ssc,'linear',2);

adcp(kk).sed_flux = fillmissing(adcp(kk).sed_flux,'linear',2);
adcp(kk).sal_flux = fillmissing(adcp(kk).sal_flux,'linear',2);
adcp(kk).sed_flux = fillmissing(adcp(kk).sed_flux,'nearest',1);
adcp(kk).sal_flux = fillmissing(adcp(kk).sal_flux,'nearest',1);

% nan values below the riverbed
for jj=1:length(adcp(kk).time)
    adcp(kk).ssc(intrp_depth>adcp(kk).depth(jj),jj) = NaN;
    adcp(kk).sal(intrp_depth>adcp(kk).depth(jj),jj) = NaN;
    adcp(kk).sed_flux(adcp(kk).z(:,1)>adcp(kk).depth(jj),jj) = NaN;
    adcp(kk).sal_flux(adcp(kk).z(:,1)>adcp(kk).depth(jj),jj) = NaN;
end

binsize = 0.5*1; %(50cm high and 1 m width assumption)% unit of kg/m3 * m/s = kg/m2/s
adcp(kk).sed_flux_east = binsize.*adcp(kk).sed_flux.*adcp(kk).east./100; 
adcp(kk).sed_flux_north = binsize.*adcp(kk).sed_flux.*adcp(kk).north./100;
adcp(kk).sal_flux_east = adcp(kk).sal_flux.*adcp(kk).east;
adcp(kk).sal_flux_north = adcp(kk).sal_flux.*adcp(kk).north;

% figure;pcolor(adcp(kk).sal),shading flat,colorbar
end



save('mmi_discharge\br_int_march2018.mat','adcp')

%% pull out a trace of the large channel for general usage:
clear all,close all,clc

cd C:\GLOVER\output\myanmar
load('mmi_discharge\br_int_march2018.mat')


ch.lat = adcp(3).lat;
ch.lon = adcp(3).lon;
ch.lat(13428:14344)=NaN;
ch.lon(13428:14344)=NaN;
[ch.x(1),ch.y(1),~] = deg2utm(ch.lat(1), ch.lon(1));
ch.dist(1)=0;
for jj=2:length(ch.lat)
    [ch.x(jj),ch.y(jj),~] = deg2utm(ch.lat(jj), ch.lon(jj));
    ch.dist(jj) = sqrt(abs(ch.x(jj)-ch.x(jj-1)).^2 + abs(ch.y(jj)-ch.y(jj-1)).^2);
end
ch.dist=cumsum(ch.dist,'omitnan');
mid=ch.dist(end)/2;
ch.dist_in=mid-(ch.dist-mid);
figure;scatter(ch.lon,ch.lat,[],ch.dist_in),colorbar

save('mmi_discharge\channel_trace.mat','ch')