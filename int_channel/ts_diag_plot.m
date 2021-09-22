%% ts_diag_plot

clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('mmi_discharge\ctdprof_8Mar18.mat');lf=ctd;
lf(11).Salinity([24:25])=NaN;
load('mmi_discharge\ctdprof_13Sep17.mat');hf=ctd;
load('mmi_discharge\channel_trace.mat')

% plot ctd locs along channel:
hflat=vertcat(hf.lat);
hflong=vertcat(hf.long);
lflat=vertcat(lf.lat);
lflong=vertcat(lf.long);

% generating background density contours
% generating background density contours
smin = 0;% set min and max values for your plot for sal and temp
smax = 18.5;
thetamin = 26.5;
thetamax = 30.5;
xdim=round((smax-smin)./0.1+1);
ydim=round((thetamax-thetamin)+1);
dens=zeros(ydim,xdim);
thetai=((1:ydim)-1)*1+thetamin;
si=((1:xdim)-1)*0.1+smin;
for j=1:ydim
    for i=1:xdim
        dens(j,i)=SW_Density(thetai(j),'C',si(i),'ppt',0.1,'bar')-1000;
    end
end


figure;
subplot(2,2,1)
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
for jj=1:length(hf)
    hftvec(jj) = 24*60*(hf(jj).time(1)-hf(1).time(1));
    d_in = hf(jj).dist_in*ones(size(hf(jj).Salinity));
    scatter(hf(jj).Salinity,hf(jj).Temperature,[],d_in,'d','filled')
end
for jj=1:length(lf)
    lftvec(jj) = 24*(lf(jj).time(1)-lf(1).time(1));
    d_in = lf(jj).dist_in*ones(size(lf(jj).Salinity));
    scatter(lf(jj).Salinity,lf(jj).Temperature,[],d_in,'*')
end
colorbar,caxis([0 6000])
colormap(cmocean('thermal'))
% ax=gca;ax.Colormap=cmocean('thermal');
xlabel('Salinity'),ylabel('Temp (^oC)')

hfdist=vertcat(hf.dist_in);
lfdist=vertcat(lf.dist_in);
subplot(222)
plot(ch.lon,ch.lat,'k-'),hold on
scatter(hflong,hflat,[],hfdist,'d','filled')
scatter(lflong,lflat,[],lfdist,'*')
caxis([0 6000])
xlabel('Lon'),ylabel('Lat')

subplot(2,2,3)
for jj=1:length(hf)
    d_in = hf(jj).dist_in*ones(size(hf(jj).Salinity));
    scatter(hf(jj).Density,hf(jj).Depth,[],d_in,'d','filled'),hold on
end
for jj=1:length(lf)
    d_in = lf(jj).dist_in*ones(size(lf(jj).Salinity));
    scatter(lf(jj).Density,lf(jj).Depth,[],d_in,'*'),hold on
end
caxis([0 6000]),axis ij
colormap(cmocean('thermal'))
xlabel('depth (m)'),ylabel('Salinity')

subplot(2,2,4)
for jj=1:length(hf)
    d_in = hf(jj).dist_in*ones(size(hf(jj).Salinity));
    scatter(hf(jj).SSCCal,hf(jj).Depth,[],d_in,'d','filled'),hold on
end
for jj=1:length(lf)
    d_in = lf(jj).dist_in*ones(size(lf(jj).Salinity));
    scatter(lf(jj).SSCCal,lf(jj).Depth,[],d_in,'*'),hold on
end
caxis([0 6000]),axis ij
xlim([0 850])
colormap(cmocean('thermal'))
xlabel('depth (m)'),ylabel('ssc (mg/L)')



%%
close all
figure;
% subplot(121)
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
for jj=1:length(hf)
    sal = nanmean(hf(jj).Salinity);
    temp = nanmean(hf(jj).Temperature);
    ssc = nanmean(lf(jj).SSCCal);
    scatter(sal,temp,ssc,hf(jj).dist/1000,'d','filled')
    scatter(sal,temp,ssc,hf(jj).dist/1000,'kd')
end
for jj=1:length(lf)
    sal = nanmean(lf(jj).Salinity);
    temp = nanmean(lf(jj).Temperature);
    ssc = nanmean(lf(jj).SSCCal);
        scatter(sal,temp,ssc,lf(jj).dist/1000,'o','filled')
        scatter(sal,temp,ssc,lf(jj).dist/1000,'ko')
end
colorbar,caxis([0 14])
C =cmocean('curl');
C = C(15:end-15,:);
colormap(C)
xlabel('Salinity'),ylabel('Temperature (^oC)')
title('Depth Avg CTD Profiles')
% 
% hfdist=vertcat(hf.dist_in);
% lfdist=vertcat(lf.dist_in);
% subplot(122)
% plot(ch.lon,ch.lat,'k-'),hold on
% scatter(hflong,hflat,[],hfdist,'d','filled')
% scatter(lflong,lflat,[],lfdist,'*')
% caxis([0 6000])
