%% TS diagrams:
% first look at just internal island ctd casts and aqd casts
% then look at ctds in river

clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('mmi_discharge\ctdprof_13Sep17.mat')
load('aqd\aqd_sep17_hc.mat'),hc=aqd; hc.sal = zeros(size(hc.time));
load('aqd\aqd_sep17_lc.mat'),lc=aqd;clear aqd
load('longterminst\BogaleRiverInstruments.mat')
load('longterminst\BogaleRiverWeather.mat')
load('mmi_discharge\channel_trace.mat')

hc.rho = SW_Density(hc.temp,'C',hc.sal,'ppt',hc.pres./10,'bar');
lc.rho = SW_Density(lc.temp,'C',lc.sal,'ppt',lc.pres./10,'bar');

% generating background density contours
smin = 0;% set min and max values for your plot for sal and temp
smax = 5;
thetamin = 26.5;
thetamax = 30;
xdim=round((smax-smin)./0.1+1);
ydim=round((thetamax-thetamin)+1);
dens=zeros(ydim,xdim);
thetai=((1:ydim)-1)*1+thetamin;
si=((1:xdim)-1)*0.1+smin;
for j=1:ydim
    for i=1:xdim
        dens(j,i)=SW_Density(thetai(j),'C',si(i),'ppt',0.1,'bar');
    end
end
dens=dens-1000;

%%
figure;

% row 1 = slope
% row 2 = ssc
% col 1 = HC LF
% col 2 = LC LF
% col 3 = LC HF

idx2 = 2740:2740+60*24;
% idx2 = 1:1+60*24;
subplot(233)% TS diagrams, colored by slope
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(lc.sal(idx2),lc.temp(idx2),[],lc.slope(idx2),'.')
title('LC hf')

subplot(236)% TS diagrams, colored by ssc
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(lc.sal(idx2),lc.temp(idx2),[],lc.ssc1(idx2),'.')
title('LC hf')

figure(2),plot(lc.time(idx2),lc.depth(idx2)),datetick('x','HH:MM','keeplimits')
%% TS diagrams: March 2018 low flow
% first look at just internal island ctd casts and aqd casts
% then look at ctds in river

clear all
cd C:\GLOVER\output\myanmar
load('mmi_discharge\ctdprof_8Mar18.mat')
load('aqd\aqd_mar18_hc.mat'),hc=aqd;
load('aqd\aqd_mar18_lc.mat'),lc=aqd;clear aqd
load('longterminst\BogaleRiverInstruments.mat')
load('longterminst\BogaleRiverWeather.mat')
load('mmi_discharge\channel_trace.mat')
ctd(11).Salinity([24:25])=NaN;
hc.rho = SW_Density(hc.temp,'C',hc.sal,'ppt',hc.pres./10,'bar');
lc.rho = SW_Density(lc.temp,'C',lc.sal,'ppt',lc.pres./10,'bar');

% generating background density contours
smin = 13.5;% set min and max values for your plot for sal and temp
smax = 18.5;
thetamin = 26.5;
thetamax = 30;
xdim=round((smax-smin)./0.1+1);
ydim=round((thetamax-thetamin)+1);
dens=zeros(ydim,xdim);
thetai=((1:ydim)-1)*1+thetamin;
si=((1:xdim)-1)*0.1+smin;
for j=1:ydim
    for i=1:xdim
        dens(j,i)=SW_Density(thetai(j),'C',si(i),'ppt',0.1,'bar');
    end
end
dens=dens-1000;

idx1 = 1:294;
% idx2 = 1:294
idx2 = find(lc.time>hc.time(1) & lc.time<hc.time(295));
tmin = [(0:5:5*357)/60]';
figure(1)
subplot(231)% TS diagrams, colored by slope
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(hc.sal(idx1),hc.temp(idx1),[],hc.slope(idx1),'.')
title('HC lf')

subplot(234)% TS diagrams, colored by slope
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(hc.sal(idx1),hc.temp(idx1),[],hc.ssc1(idx1),'.')
title('HC lf')

subplot(232)% TS diagrams, colored by slope
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(lc.sal(idx2),lc.temp(idx2),[],lc.slope(idx2),'.')
title('LC lf')

subplot(235)% TS diagrams, colored by ssc
[c,h]=contour(si,thetai,dens,'k');
hold on
clabel(c,h,'LabelSpacing',1000);
scatter(lc.sal(idx2),lc.temp(idx2),[],lc.ssc1(idx2),'.')
title('LC lf')



figure(3),plot(lc.time(idx2),lc.depth(idx2)),datetick('x','HH:MM','keeplimits')
figure(4),plot(hc.time(idx1),hc.depth(idx1)),datetick('x','HH:MM','keeplimits')
%%
figure(1)
for jj=1:3
    subplot(2,3,jj)
    colorbar,caxis([-0.8 0.8])
    xlabel('Salinity'),ylabel('Temp (^oC)')
    s=gca; s.Colormap=cmocean('balance');
    ylim([26.5 30.5])
end
 
for jj=4:6
    subplot(2,3,jj)
    colorbar,caxis([0 150])
    xlabel('Salinity'),ylabel('Temp (^oC)')
    s=gca; s.Colormap=cmocean('amp');
    ylim([26.5 30.5])
end
      
