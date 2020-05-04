

clear all%,close all,clc
load('AyeMar18_CTD_all.mat')
ctd=CombinedProfiles; clear CombinedProfiles

% Spring = 5 March 2018; Neap = 9 March 2018

% find distance along transect
% find idealized transect end points (from adcp processing:
% BR_neap_Residuals)
[right.x,right.y]=deg2utm(16.099239, 95.320838);
[left.x,left.y]=deg2utm(16.100713, 95.330033);
cc=1;
for jj=1:length(ctd)   
    % convert the coordinates of your ADCP records to UTM
    [ctd(jj).x, ctd(jj).y,~] = deg2utm(ctd(jj).lat, ctd(jj).long);
    % project the ADCP data onto the idealized transect
    ctdproj = proj([left.x - ctd(jj).x, left.y - ctd(jj).y],...
        [left.x- right.x, left.y - right.y]);
    % define new distances along the idealized transect
    ctd(jj).dist =  round(sqrt(ctdproj(1).^2 + ctdproj(2).^2));
    if ctd(jj).dist>300 && ctd(jj).dist<650
        trans(cc)=ctd(jj);
        cc=cc+1;
    end
end

for jj=1:length(trans)
    tvec(jj)=trans(jj).time(1);
end
[tvec,idx]=sort(tvec,'ascend');
tvec=datevec(tvec);tvec(:,6)=0;
daychk=tvec(:,3);
tvec=datenum(tvec);
trans=trans(idx);


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

for kk=[5,9]
dd = trans(daychk==kk);
tvec_dd=tvec(daychk==kk,:);

L=round((tvec_dd(end)-tvec_dd(1))*60*24+1);
sigmaDepth=linspace(0,1,1000);
sigmaSalinity=NaN([1000,L]);
sigmaSSC=NaN([1000,L]);

for jj=1:length(dd)
    tidx=round((tvec_dd(jj)-tvec_dd(1))*60*24+1);
    dd(jj).Depth=dd(jj).Depth./max(dd(jj).Depth);
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    sigmaSalinity(:,tidx)=interp1(dd(jj).Depth,dd(jj).Salinity,sigmaDepth);
    sigmaSSC(:,tidx)=interp1(dd(jj).Depth,dd(jj).SSCCal,sigmaDepth);
end
% interpolate to surface
sigmaSalinity=fillmissing(sigmaSalinity,'nearest',1);
sigmaSSC=fillmissing(sigmaSSC,'nearest',1);
% interpolate horizontally between casts
sigmaSalinity=fillmissing(sigmaSalinity,'linear',2);
sigmaSSC=fillmissing(sigmaSSC,'linear',2);

timeplot=tvec_dd(1):(1/(24*60)):tvec_dd(end);

% plot the salinity
levels=0:1:20;
figure;
subplot(2,1,1)
contourf(timeplot,sigmaDepth,sigmaSalinity,levels,'showtext','on')
colormap(C1),shading flat,axis ij,cbar=colorbar;caxis([0 11])
ylabel(cbar,'Salinity')
ylabel('Depth')
datetick('x','HHMM','keeplimits')

% plot the ssc
levels=0:100:400;
subplot(2,1,2)
contourf(timeplot,sigmaDepth,sigmaSSC,levels,'showtext','on')
% pcolor(1:60,sigmaDepth,sigmaSSC)
shading flat,axis ij,cbar=colorbar;caxis([0 400])
ylabel(cbar,'mg/L'),ylabel('Depth')
datetick('x','HHMM','keeplimits'),xlabel(datestr(round(tvec_dd(1))))

end
