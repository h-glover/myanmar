
clear all%,close all,clc
load('AyeMar18_CTD_all.mat')
for jj=1:length(CombinedProfiles)
    idx(jj)=CombinedProfiles(jj).time(1);
end
[time,idx]=sort(idx,'ascend');
ctd=CombinedProfiles(idx);

for jj=1:length(ctd)
    ctd(jj).profileidx=jj; %label 1:length(ctd) for later use
    ctd(jj).timeidx=ctd(jj).time(1);%time cast starts, for matching adcp transect
end
% concatenate the profile indexs and times
allidx_ctd=horzcat(ctd.profileidx);
alltime_ctd=vertcat(ctd.timeidx);
alltime_ctd=datevec(alltime_ctd);alltime_ctd(:,6)=0;
alltime_ctd=datenum(alltime_ctd);

load('BR_Sept17_Spring.mat')
% adcp time and index will be reused in the SedFlux calculation
for jj=1:length(adcp)
    %get the matlab time for each transect and assign a transect number to
    %each time point
    adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
        adcp(jj).hour,adcp(jj).minute,adcp(jj).second)+(6.5/24);
    adcp(jj).transidx=jj*ones(1,length(adcp(jj).time));
    idx(jj)=adcp(jj).time(1);
end

% concatenate into 1 vector 
allidx_adcp=horzcat(adcp.transidx);
alltime_adcp=horzcat(adcp.time);
alltime_adcp=datevec(alltime_adcp);alltime_adcp(:,6)=0;
alltime_adcp=datenum(alltime_adcp);

% compare adcp time and ctd time (seconds removed)
[~,iA,iC]=intersect(alltime_adcp,alltime_ctd);

figure;
plot(alltime_ctd,ones(1,length(alltime_ctd)),'ok'),hold on
plot(alltime_ctd(iC),ones(1,length(iC)),'r*')
plot(alltime_adcp,allidx_adcp-2,'dk')
% % pull out all ctd profiles from this day
% ctd_trans=allidx_ctd(iC);
% allidx_adcp=allidx_adcp(iA);
% ctd=ctd(iC);