clear all,close all,clc
% % define colormap:
% C1=linspace(0.6,0.95,400)';
% C1(:,2)=linspace(0.1,0.95,400)';
% C1(:,3)=linspace(0.1,1,400)';
% C1=flipud(C1);
% C2=linspace(0.4,1,500)';
% C2(:,2)=linspace(0.3,0.9,500)';
% C2(:,3)=linspace(0.5,0.9,500)';
% C2=flipud(C2);

load('AyeMar18_CTD_all.mat')
for jj=1:length(CombinedProfiles)
    daychk(jj,:)=datevec(CombinedProfiles(jj).time(1));
end
daychk=daychk(:,3);
[daychk,idx]=sort(daychk);
CombinedProfiles=CombinedProfiles(idx);
CombinedProfiles([9,12,17])=[];
daychk([9,12,17])=[];% remove repeat profiles
figure;
for kk=1:2 %do for both days on Yangon

% dd = CombinedProfiles(daychk==1 | daychk==2);
dd = CombinedProfiles(daychk==kk);

% project the ctd locations to the thalweg to get distance along estuary
lat=vertcat(dd.lat);
lon=vertcat(dd.long);
[loc.x,loc.y]=deg2utm(lat,lon);

% define 3 transects along estuary for distance along river calculation
[mouth.x,mouth.y]=deg2utm(16.447348, 96.352243);
[split.x,split.y]=deg2utm( 16.756927, 96.199459);
[turn.x,turn.y]=deg2utm( 16.775776, 96.119558);
[top.x,top.y]=deg2utm(16.880264,96.087471);
% project the casts onto the nearest transect along the river
for jj=1:length(dd)
    if loc.x(jj)>=split.x %lower river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);
    elseif loc.x(jj)>=turn.x && loc.x(jj)<split.x %middle river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);    
    elseif loc.x(jj)<turn.x %upper river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);    
    end
end
% define new distances along the idealized transect 0=mouth
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2)/1000);

% convert depth to sigma coordiantes and put the salinity profiles in the
% right locations along the river
sigmaDepth=linspace(0,1,1000);
sigmaSalinity=NaN([length(sigmaDepth),60]);
sigmaSSC=NaN([length(sigmaDepth),60]);
for jj=1:length(dd)
    dd(jj).Depth=dd(jj).Depth./max(dd(jj).Depth);
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    sigmaSalinity(:,dist(jj))=interp1(dd(jj).Depth,dd(jj).Salinity,sigmaDepth);
    sigmaSSC(:,dist(jj))=interp1(dd(jj).Depth,dd(jj).SSCCal,sigmaDepth);

end
% interpolate to surface
sigmaSalinity=fillmissing(sigmaSalinity,'nearest',1);
sigmaSSC=fillmissing(sigmaSSC,'nearest',1);
% interpolate horizontally between casts
% interpolate horizontally between casts
fillidx=min(dist):max(dist);
sigmaSalinity(:,fillidx)=fillmissing(sigmaSalinity(:,fillidx),'linear',2);
sigmaSSC(:,fillidx)=fillmissing(sigmaSSC(:,fillidx),'linear',2);

% plot the salinity
levels=0:1:14;
figure(1)
subplot(2,1,kk)
contourf(1:60,sigmaDepth,sigmaSalinity,levels,'showtext','on')
colormap(C1),shading flat,axis ij,cbar=colorbar;caxis([0 14])
ylabel(cbar,'Salinity')
ylabel('Depth')
title(['Salinity (upriver, against spring ebb) ',datestr(dd(1).time(1),'mmm-dd-yyyy, HHMM')])


% plot the ssc
levels=0:500:4000;
figure(2)
subplot(2,1,kk)
contourf(1:60,sigmaDepth,sigmaSSC,levels,'showtext','on')
% pcolor(1:60,sigmaDepth,sigmaSSC)
shading flat,axis ij,cbar=colorbar;colormap(C2),caxis([0 3500])
ylabel('Depth'),ylabel(cbar,'mg/L')
title(['SSC (upriver, against spring ebb) ',datestr(dd(1).time(1),'mmm-dd-yyyy, HHMM')])

end
xlabel('km from mouth')
figure(1),xlabel('km from mouth')


%% VIMS profiles Dec 2017

clear all%,close all,clc
% define colormap:
C1=linspace(0.6,0.95,400)';
C1(:,2)=linspace(0.1,0.95,400)';
C1(:,3)=linspace(0.1,1,400)';
C1=flipud(C1);
C2=linspace(0.4,1,500)';
C2(:,2)=linspace(0.3,0.9,500)';
C2(:,3)=linspace(0.5,0.9,500)';
C2=flipud(C2);

load('AyeDec17_CTD_all.mat')


% dd = CombinedProfiles(daychk==1 | daychk==2);
dd = profiles;

% project the ctd locations to the thalweg to get distance along estuary
lat=vertcat(dd.lat);
lon=vertcat(dd.long);
[loc.x,loc.y]=deg2utm(lat,lon);

% define 3 transects along estuary for distance along river calculation
[mouth.x,mouth.y]=deg2utm(16.447348, 96.352243);
[split.x,split.y]=deg2utm( 16.756927, 96.199459);
[turn.x,turn.y]=deg2utm( 16.775776, 96.119558);
[top.x,top.y]=deg2utm(16.880264,96.087471);
% project the casts onto the nearest transect along the river
for jj=1:length(dd)
    if loc.x(jj)>=split.x %lower river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);
    elseif loc.x(jj)>=turn.x && loc.x(jj)<split.x %middle river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);    
    elseif loc.x(jj)<turn.x %upper river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);    
    end
end
% define new distances along the idealized transect 0=mouth
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2)/1000);

% convert depth to sigma coordiantes and put the salinity profiles in the
% right locations along the river
sigmaDepth=linspace(0,1,1000);
sigmaSalinity=NaN([length(sigmaDepth),60]);
sigmaSSC=NaN([length(sigmaDepth),60]);
for jj=1:length(dd)
    dd(jj).Depth=dd(jj).Depth./max(dd(jj).Depth);
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    sigmaSalinity(:,dist(jj))=interp1(dd(jj).Depth,dd(jj).Salinity,sigmaDepth);
    sigmaSSC(:,dist(jj))=interp1(dd(jj).Depth,dd(jj).SSC,sigmaDepth);

end
% interpolate to surface
sigmaSalinity=fillmissing(sigmaSalinity,'nearest',1);
sigmaSSC=fillmissing(sigmaSSC,'nearest',1);
% interpolate horizontally between casts
% interpolate horizontally between casts
fillidx=min(dist):max(dist);
sigmaSalinity(:,fillidx)=fillmissing(sigmaSalinity(:,fillidx),'linear',2);
sigmaSSC(:,fillidx)=fillmissing(sigmaSSC(:,fillidx),'linear',2);

% plot the salinity
levels=0:1:16;
figure;
contourf(1:60,sigmaDepth,sigmaSalinity,levels,'showtext','on')
colormap(C1),shading flat,axis ij,colorbar,caxis([0 16])
ylabel('Depth')
title(['Salinity (upriver, with flood) ',datestr(dd(1).time(1),'mmm-dd-yyyy, HHMM')])


% % plot the ssc
% levels=0:500:4000;
% figure;
% contourf(1:60,sigmaDepth,sigmaSSC,levels,'showtext','on')
% % pcolor(1:60,sigmaDepth,sigmaSSC)
% shading flat,axis ij,colorbar,colormap(C2),caxis([0 3500])
% ylabel('Depth')
% title(['SSC (upriver, with flood) ',datestr(dd(1).time(1),'mmm-dd-yyyy, HHMM')])
% 

xlabel('km from mouth')
figure(1),xlabel('km from mouth')

% %% 
% load('YangonRiverInstruments.mat')
% idx=find(YangonRiver.datenum>profiles(1).time(1)-5/24 &...
%     YangonRiver.datenum<profiles(end).time(end)+5/24);
% figure;
% plot(YangonRiver.datenum(idx),YangonRiver.Depth(idx))
% 
% 
% % VIMS profile locations
% figure;
% plot(lon,lat,'k*')
% hold on
% plot(96.352243,16.447348,'rd')%mouth
% plot(96.199459,16.756927,'r^')%Yangon split

%%
clear all,close all,clc
load('AyeSept17_CTD_all.mat')

% find distance along transect
% find idealized transect end points (from adcp processing:
% BR_neap_Residuals)
profiles=CombinedProfiles;

for jj=1:length(profiles)
    tvec(jj)=profiles(jj).time(1);
end
[tvec,idx]=sort(tvec,'ascend');
tvec=datevec(tvec);tvec(:,6)=0;
daychk=tvec(:,3);
tvec=datenum(tvec);
profiles=profiles(idx);


%
% convert depth to sigma coordiantes and put the salinity profiles in the
% right locations along the river
C1=linspace(0.6,0.95,400)';
C1(:,2)=linspace(0.1,0.95,400)';
C1(:,3)=linspace(0.1,1,400)';
C1=flipud(C1);
C2=linspace(0.4,1,500)';
C2(:,2)=linspace(0.3,0.9,500)';
C2(:,3)=linspace(0.5,0.9,500)';
C2=flipud(C2);

kk=17;
dd = profiles(daychk==kk);
tvec_dd=tvec(daychk==kk,:);

% project the ctd locations to the thalweg to get distance along estuary
lat=vertcat(dd.lat);
lon=vertcat(dd.long);
[loc.x,loc.y]=deg2utm(lat,lon);

% define 3 transects along estuary for distance along river calculation
[mouth.x,mouth.y]=deg2utm(16.447348, 96.352243);
[split.x,split.y]=deg2utm( 16.756927, 96.199459);
[turn.x,turn.y]=deg2utm( 16.775776, 96.119558);
[top.x,top.y]=deg2utm(16.880264,96.087471);
% project the casts onto the nearest transect along the river
for jj=1:length(dd)
    if loc.x(jj)>=split.x %lower river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);
    elseif loc.x(jj)>=turn.x && loc.x(jj)<split.x %middle river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);    
    elseif loc.x(jj)<turn.x %upper river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);    
    end
end
% define new distances along the idealized transect 0=mouth
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2)/1000);

L=round((tvec_dd(end)-tvec_dd(1))*60*24+1);
sigmaDepth=linspace(0,1,100);
sigmaSSC=NaN([length(sigmaDepth),60]);
for jj=1:length(dd)
    dd(jj).Depth=dd(jj).Depth./max(dd(jj).Depth);
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    sigmaSSC(:,dist(jj))=interp1(dd(jj).Depth,dd(jj).SSCCal,sigmaDepth);

end
% interpolate to surface
sigmaSSC=fillmissing(sigmaSSC,'nearest',1);
% interpolate horizontally between casts
fillidx=min(dist):max(dist);
sigmaSSC(:,fillidx)=fillmissing(sigmaSSC(:,fillidx),'linear',2);


figure;
% plot the ssc
levels=0:100:600;
contourf(1:60,sigmaDepth,sigmaSSC,levels,'showtext','on')
colormap(C2)
shading flat,axis ij,cbar=colorbar;caxis([0 600])
ylabel(cbar,'mg/L'),ylabel('Depth')
title('Yangon, High Flow, Neap')
xlabel('km from mouth'),xlim([15 45])