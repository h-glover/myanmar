
%% ALL THIS IS NOW IN MMI_DISCHARGE_CALC CODE!!!
% calc richardson number at each ctd station by pulling out the velocity
% from the adcp

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

% calculate density and density gradient for each depth:
% figure;
for jj=1:length(ctd)
    ctd(jj).Density = sw_dens(ctd(jj).Salinity,ctd(jj).Temperature,ctd(jj).Depth);
    % fix Salinity blips at bottom and top of some profiles
    if jj>14 || jj<5
        ctd(jj).Density(end-2:end) = ctd(jj).Density(end-4);
        ctd(jj).Density(1:2) = ctd(jj).Density(3);
    end
    ctd(jj).Density = fillmissing(ctd(jj).Density,'nearest');
%     ctd(jj).Density = movmean(ctd(jj).Density,3);
    
    % calc depths and density gradient for in between each point:
    ctd(jj).dens_grad = interp1(ctd(jj).Depth,ctd(jj).Density,adcp(3).z(:,1));
    
    
    ctd(jj).N2 = (ctd(jj).dens_grad(1:end-2) - ctd(jj).dens_grad(3:end))./1;
    ctd(jj).N2 = [ctd(jj).N2(1);ctd(jj).N2;ctd(jj).N2(end)];
    ctd(jj).N2 = -9.8/1010*ctd(jj).N2;
    ctd(jj).N2 = abs(ctd(jj).N2);
 
%     subplot(4,5,jj),plot(ctd(jj).N2,adcp(3).z(:,1))
    
    % get the adcp transect
    ctd(jj).transect=allidx_adcp(jj);
end


% map ctd (N2 vals) onto adcp transect3
adcp(3).N2 = NaN(size(adcp(3).spd));
[~,iC,iA] = intersect(alltime_ctd,adcp(3).time);
for jj=1:length(iC)
    adcp(3).N2(:,iA(jj)) = ctd(jj).N2;
end
adcp(3).N2 = fillmissing(adcp(3).N2,'nearest',2);

adcp(3).KE = [adcp(3).spd(1,:);adcp(3).spd;adcp(3).spd(end,:)]./100;
adcp(3).KE = movmean(adcp(3).KE,5);
adcp(3).KE = movmean(adcp(3).KE,7,2);
adcp(3).KE = (adcp(3).KE(1:end-2,:) - adcp(3).KE(3:end,:)).^2;
adcp(3).KE(adcp(3).KE>0.4) = NaN;

adcp(3).Ri = adcp(3).N2./adcp(3).KE;


% adcp(3).Ri = fillmissing(adcp(3).Ri,'nearest',2);
% adcp(3).Ri = movmean(adcp(3).Ri,3,2);
% adcp(3).Ri = movmean(adcp(3).Ri,3);
% 
% for jj=1:length(adcp(3).time)
%     adcp(3).Ri(adcp(3).z(:,1)>adcp(3).depth(jj),jj) = NaN;
% end   
    
    
figure;
subplot(411)
pcolor(adcp(3).time,adcp(3).z(:,1),adcp(3).spd),caxis([0 100])
subplot(412)
pcolor(adcp(3).time,adcp(3).z(:,1),adcp(3).N2),caxis([0 0.01])
subplot(413)
pcolor(adcp(3).time,adcp(3).z(:,1),adcp(3).KE),caxis([0 0.01])
subplot(414)
pcolor(adcp(3).time,adcp(3).z(:,1),adcp(3).Ri),caxis([0 2]),hold on
ax = gca; ax.Colormap = cmocean('balance','pivot',0.25);
% levels=[0,0.25];
% hold on
% [t,ia,ic] = unique(adcp(3).time);
% contour(t,adcp(3).z(:,1),adcp(3).Ri(:,ia),levels,'Color','k')

for jj=1:4
    subplot(4,1,jj),shading flat,colorbar
    axis ij
    ylim([0 15])
end

