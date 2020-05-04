
clear all,close all,clc
% define colormap:
C1=linspace(0.6,0.95,400)';
C1(:,2)=linspace(0.1,0.95,400)';
C1(:,3)=linspace(0.1,1,400)';
C1=flipud(C1);
C2=linspace(0.4,1,500)';
C2(:,2)=linspace(0.3,0.9,500)';
C2(:,3)=linspace(0.5,0.9,500)';
C2=flipud(C2);

load('AyeMar18_CTD_all.mat')
for jj=1:length(CombinedProfiles)
    daychk(jj)=CombinedProfiles(jj).time(1);
end
[daychk,idx]=sort(daychk,'ascend');
daychk=datevec(daychk);
daychk=daychk(:,3);
CombinedProfiles=CombinedProfiles(idx);


cc=1;
for kk=[14,15]
dd = CombinedProfiles(daychk==kk);


% project the ctd locations to the thalweg to get distance along estuary
lat=vertcat(dd.lat);
lon=vertcat(dd.long);
[loc.x,loc.y]=deg2utm(lat,lon);

% define 3 transects along estuary for distance along river calculation
[mouth.x,mouth.y]=deg2utm(16.354683, 94.676328);
[top.x,top.y]=deg2utm(16.778223,94.728839);
% project the casts onto the nearest transect along the river
for jj=1:length(dd)
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- top.x, mouth.y- top.y]);
end
% define new distances along the idealized transect 0=mouth
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2)/1000)+58;

% convert depth to sigma coordiantes and put the salinity profiles in the
% right locations along the river
sigmaDepth=linspace(0,1,1000);
sigmaSalinity=NaN([length(sigmaDepth),100]);
sigmaSSC=NaN([length(sigmaDepth),100]);
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
fillidx=min(dist):max(dist);
sigmaSalinity(:,fillidx)=fillmissing(sigmaSalinity(:,fillidx),'linear',2);
sigmaSSC(:,fillidx)=fillmissing(sigmaSSC(:,fillidx),'linear',2);

% plot the salinity
levels=0:1:20;
figure(1)
subplot(2,1,cc)
contourf(1:100,sigmaDepth,sigmaSalinity,levels,'showtext','on')
colormap(C1),shading flat,axis ij,cbar=colorbar;caxis([0 12])
ylabel(cbar,'Salinity')
ylabel('Depth'),xlim([50 100])
title(['Salinity (downriver against flood) ',datestr(dd(1).time(1),'mmm-dd-yyyy, HHMM')])


% plot the ssc
levels=0:250:4000;
figure(2)
subplot(2,1,cc)
contourf(1:100,sigmaDepth,sigmaSSC,levels,'showtext','off')
% pcolor(1:60,sigmaDepth,sigmaSSC)
shading flat,axis ij,cbar=colorbar;colormap(C2),caxis([0 500])
ylabel('Depth'),xlim([50 100]),ylabel(cbar,'mg/L')
title(['SSC (downriver against flood) ',datestr(dd(1).time(1),'mmm-dd-yyyy, HHMM')])
cc=1+cc;
end
xlabel('km from mouth')
figure(1),xlabel('km from mouth')
