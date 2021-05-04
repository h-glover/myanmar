% calc inundation period for croc and cyclone
clear all,close all,clc
cd C:\GLOVER\output\myanmar\
% compare elevation to water level
load('SurveyPts2019.mat')
load('longterminst\BogaleRiverInstruments.mat')
br = BogaleRiver; clear BogaleRiver

% sept2019 aquadopp data
load('aqd\aqd_sep19_hc.mat'),hc=aqd;
load('aqd\aqd_sep19_lc.mat'),lc=aqd;
load('aqd\aqd_sep19_ag.mat'),ag=aqd;

% fix depth in high conn:
hc.depth_elev(387:end) = hc.depth_elev(387:end) + 0.3;
hc.depth_elev(1108:end) = hc.depth_elev(1108:end) + 0.05;

%% boxplot of land surface elevations
boxdat.Z = survey.Zcorr_wl(survey.ID==1 | survey.ID==2 | survey.ID==6 & survey.type==0);
boxdat.Z(boxdat.Z<0 | boxdat.Z>3)=NaN;
boxdat.ID = survey.ID(survey.ID==1 | survey.ID==2 | survey.ID==6 & survey.type==0);

figure;
boxplot(boxdat.Z,boxdat.ID,'GroupOrder',{'2','1','6'},'Symbol','k')
ax=gca; ax.XTickLabel = {'Cyclone','Croc Stn','Agri Field'};
ylabel('Elevation (m)')
%% calc inundation period for croc

figure;
plot(br.datenum,br.MeinmahlaWaterLevel,'b'),hold on
plot(lc.time,lc.depth_elev,'b')
R=refline(0,lc.flat); R.Color='b';
plot(ag.time,ag.depth_elev,'k')
R=refline(0,ag.flat); R.Color='k';
plot(hc.time,hc.depth_elev,'r')
R=refline(0,hc.flat); R.Color='r';
xlim([ag.time(1)-5 ag.time(end)])
datetick('x','dd','keeplimits')
% ag.flat = nanmean(survey.Zcorr_wl([3:15,34:54]));
% hc.flat = nanmean(survey.Zcorr_wl(161:176));
% lc.flat= nanmean(survey.Zcorr_wl(126:139));
% 
% ag.flat - max(ag.depth_elev)
% hc.flat - max(hc.depth_elev)

idx = find(ag.depth_elev(500:end)>ag.flat);
t_inun = (length(idx)*(ag.time(2)-ag.time(1))*24*60)/4;
idx = find(hc.depth_elev>hc.flat);
t_inun = (length(idx)*(hc.time(2)-hc.time(1))*24*60)/4;
idx = find(lc.depth_elev>lc.flat);
t_inun = (length(idx)*(lc.time(2)-lc.time(1))*24*60)/4;

