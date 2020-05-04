

clear all,close all,clc
load('AyeJan19_CTD_all.mat')
load('YangonWL_Prediction.mat')
% Spring = 21 March 2018; Neap = 14 March 2018

% find distance along transect
% find idealized transect end points (from adcp processing:
% BR_neap_Residuals)


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

for kk=[14,21]
dd = profiles(daychk==kk);
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
    sigmaSSC(:,tidx)=interp1(dd(jj).Depth,dd(jj).SSC,sigmaDepth);
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
subplot(3,1,1)
contourf(timeplot,sigmaDepth,sigmaSalinity,levels,'showtext','on')
hold on
plot(tvec_dd,zeros([1,length(tvec_dd)]),'kv')
colormap(C1),shading flat,axis ij,cbar=colorbar;caxis([0 11])
ylabel(cbar,'Salinity')
ylabel('Depth')
datetick('x','HHMM','keeplimits')

% plot the ssc
levels=0:1000:4000;
subplot(3,1,2)
contourf(timeplot,sigmaDepth,sigmaSSC,levels,'showtext','on')
hold on
plot(tvec_dd,zeros([1,length(tvec_dd)]),'kv')
shading flat,axis ij,cbar=colorbar;caxis([0 4000])
ylabel(cbar,'mg/L'),ylabel('Depth')
datetick('x','HHMM','keeplimits')

subplot(313)
plot(p_time,YOUT2)
xlim([tvec_dd(1) tvec_dd(end)])
datetick('x','HHMM','keeplimits'),xlabel(datestr(round(tvec_dd(1))))
ylabel('water level')

end
subplot(311),title('Yangon, Spring')
figure(1),subplot(311),title('Yangon, Neap')



