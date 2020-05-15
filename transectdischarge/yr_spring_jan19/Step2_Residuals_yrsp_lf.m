clear all,close all,clc
% load the mat file containng the adcp structure produced by running Dan Nowacki's adcp_ascii function
cd C:\GLOVER\output\myanmar\transectdischarge\yr_spring_jan19

% load the spring yangon data to get the depth profile:
load('C:\GLOVER\output\myanmar\transectdischarge\yr_neap_sept17\YR_Sept17_Neap.mat')
depths_hf=vertcat(adcp.interpdepths);
clear adcp
depths_hf(depths_hf<0)=NaN;
% subtract 1 meter for shift from hf to lf
depths_hf_mean=nanmean(depths_hf)-1;
% load water level data to compare to discharge:
load('C:\GLOVER\output\myanmar\longterminst\YangonRiverInstruments.mat')
tvec=days(datenum('21-Jan-2019'):minutes(1):datenum('22-Jan-2019'));
wl=interp1(YangonRiver.datenum,YangonRiver.ut_Depth,tvec);

% load the raw enx data
load('YR_spring_jan19.mat')
c=1;
UpDists=[];
DownDists=[];
AllVels=[];
transect_length=1000; % in meters

% define the right and left bank coordinates of your idealized transect in UTM
[right.x,right.y]=deg2utm(16.758914, 96.168642); %Yangon
[left.x,left.y]=deg2utm(16.766158, 96.170355); %Yangon
% determine the heading of along-channel flow (orthogonal to transect)
% from idealized transect end point coordinates
deg2rotate=90+atand(abs(left.x-right.x)/abs(left.y-right.y));
cmap = cmocean('balance');
for n=1:length(adcp)
    
    % convert the coordinates of of your ADCP records to UTM
    [adcp(n).x, adcp(n).y,~] = deg2utm(adcp(n).lat(~isnan(adcp(n).lat)), adcp(n).lon(~isnan(adcp(n).lon)));
    
    
    % project the ADCP data onto the idealized transect
    for qq = 1:length(adcp(n).x)
        adcp(n).proj(qq,:) = proj([left.x - adcp(n).x(qq), left.y - adcp(n).y(qq)], [left.x- right.x, left.y - right.y]);
    end
    
    % define new distances along the idealized transect
    adcp(n).dist =  sqrt(adcp(n).proj(:,1).^2 + adcp(n).proj(:,2).^2);
    
    %     deg2rotate=127;
    
    % perform the coordinate rotation into along and accross components
    [adcp(n).across, adcp(n).along] = rot_earth(adcp(n).east, adcp(n).north, deg2rotate);
    
    % define transect width
    adcp(n).interpdist = [0:1:transect_length]';
    
    % find only the unique data (eliminate distance repeats)
    [adcp(n).uniquedist,uniqueindices]=unique(adcp(n).dist);
    
    % make depth vector based on mean depth
    t = datevec(adcp(n).time(1));t(6)=0; t=datenum(t);
    adcp(n).interpdepths=depths_hf_mean+wl(tvec==t);
    
    
    % interpolate the
    adcp(n).interpalong = interp1(adcp(n).uniquedist, adcp(n).along(:,uniqueindices)', adcp(n).interpdist)';
    adcp(n).interpacross = interp1(adcp(n).uniquedist, adcp(n).across(:,uniqueindices)', adcp(n).interpdist)';
    
    adcp(n).interpspd = sqrt(adcp(n).interpalong.^2 + adcp(n).interpacross.^2);
%     adcp(n).spd = sqrt(adcp(n).east.^2 + adcp(n).north.^2);
    [adcp(n).spd,adcp(n).dir]=uv2sd(adcp(n).east,adcp(n).north,adcp(n).up);
    % compute x-sectional area
%     figure;
%     subplot(3,1,1:2),pcolor(0:1000,adcp(n).bindepth,adcp(n).interpalong)
%     shading flat,colorbar,caxis([-200 200]),colormap(cmap)
%     hold on, plot(0:1000,adcp(n).interpdepths)
%     subplot(313)
%     plot(tvec,wl),hold on,plot([adcp(n).time(1),adcp(n).time(1)],[-3,3],'r-')
%     datetick('x','HH','keeplimits')
    
%     figure;
%     subplot(2,1,1),pcolor(0:1000,adcp(n).bindepth,adcp(n).spd)
%     shading flat,colorbar,caxis([-200 200]),colormap(cmap)
%     hold on, plot(0:1000,adcp(n).interpdepths)
%     subplot(212),pcolor(0:1000,adcp(n).bindepth,adcp(n).interpacross)
%     shading flat,colorbar,caxis([-200 200]),colormap(cmap)
%     hold on, plot(0:1000,adcp(n).interpdepths)

    figure;
    subplot(311),pcolor(1:length(adcp(n).time),adcp(n).bindepth,adcp(n).spd)
    shading flat,colorbar,caxis([0 300]),colormap(cmocean('tempo'))
    yyaxis right,plot(nanmean(adcp(n).spd),'r')
    ylim([0 300])
    
    subplot(312),pcolor(1:length(adcp(n).time),adcp(n).bindepth,adcp(n).dir)
    shading flat,colorbar,caxis([0 360])
    yyaxis right,plot(nanmean(adcp(n).dir),'r'),ylim([0 360])
    
    subplot(313)
    plot(tvec,wl),hold on,plot([adcp(n).time(1),adcp(n).time(1)],[-3,3],'r-')
    datetick('x','HH','keeplimits')
   
    
    c=c+1;
end

% save YR_Jan19_spring.mat adcp

%%

for jj=1:length(adcp)
    figure;
    quiver(adcp(jj).lon,adcp(jj).lat,...
        adcp(jj).east(2,:)',adcp(jj).north(2,:)')

end


%%
clear all,close all,clc
F=dir('*.mat');
% figure;
for ff=1:length(F)
load(F(ff).name)
figure;%subplot(5,1,ff)
for jj=1:4
    quiver(adcp(jj).lon,adcp(jj).lat,...
        adcp(jj).east(2,:)',adcp(jj).north(2,:)')
    hold on
end
legend
title(F(ff).name)
end
