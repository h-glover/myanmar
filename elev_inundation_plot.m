% calc inundation period for croc and cyclone
clear all,close all,clc
cd C:\GLOVER\output\myanmar\
% compare elevation to water level
load('SurveyPts2019.mat')
load('longterminst\BogaleRiverInstruments.mat')
br = BogaleRiver; clear BogaleRiver


% boxplot of land surface elevations
boxdat.Z = survey.Zcorr_wl(survey.ID==1 | survey.ID==2 | survey.ID==6 & survey.type==0);
boxdat.Z(boxdat.Z<0 | boxdat.Z>3)=NaN;
boxdat.ID = survey.ID(survey.ID==1 | survey.ID==2 | survey.ID==6 & survey.type==0);
boxdat.Z(boxdat.ID==2 & boxdat.Z>1.5)=NaN;


figure;
subplot(211)
boxplot(boxdat.Z,boxdat.ID,'GroupOrder',{'2','1','6'},'Symbol','k')
ax=gca; ax.XTickLabel = {'Cyclone','Croc Stn','Agri Field'};
ylabel('Elevation (m)'),ylim([0 2.5])




%%
clear all

cd C:\GLOVER\output\myanmar\aqd
F= dir('*.mat');

for ff=1:length(F)
    load(F(ff).name)
    tr(ff) = max(aqd.depth,[],'omitnan')-min(aqd.depth,[],'omitnan');
    mwl(ff) = nanmean(aqd.depth);
    rnge=[]; rnge = abs(aqd.depth-mwl(ff));
    avg_range(ff) = nanmean(rnge);
end
tr([7,8]) = NaN;%2.92; % manual fix for bad data in HC Sept 2019
avg_range=avg_range*2;
% distance up river, cyclone is 11
loc = [39,22,17,22,17,39,22,17,11,11];
ag = [1,6];
hc = [2,4,7];
lc = [3,5,8];

hf = [4:9];
lf = [1:3,10];


% load the cyclone shelter data:
load('C:\GLOVER\output\myanmar\longterminst\BogaleRiverInstruments.mat')
br = BogaleRiver; clear BogaleRiver
hf_dates = 132900:144100;
msl_hf = nanmean(br.MeinmahlaDepth(hf_dates));
cyc_range = abs(br.MeinmahlaDepth(hf_dates) - msl_hf);
avg_range(9) = nanmean(cyc_range)*2;
tr(9) = max(br.MeinmahlaDepth(hf_dates)) - min(br.MeinmahlaDepth(hf_dates));

lf_dates = 96421:107328;
msl_lf = nanmean(br.MeinmahlaDepth(lf_dates));
cyc_range = abs(br.MeinmahlaDepth(lf_dates) - msl_lf);
avg_range(10) = nanmean(cyc_range)*2;
tr(10) = max(br.MeinmahlaDepth(lf_dates)) - min(br.MeinmahlaDepth(lf_dates));

% plot all
subplot(212)
yyaxis left
plot(loc(hf),tr(hf)/2,'ko'),hold on
plot(loc(lf),tr(lf)/2,'ro')
% plot(loc(hf),avg_range(hf),'k*')
% plot(loc(lf),avg_range(lf),'r*')
legend({'high flow','low flow'})
% legend({'hf max range','lf max range','hf avg range','lf avg range'})
% ylim([0 3.6]),xlim([0 45])
ylabel('Max. observed tidal amp (m)'),xlabel('Distance from mouth (km)')

ylim([0 2.5])

yyaxis right
% add lines:
x = [11000,39000];
y1 = 2.2 - 0.00004*(x(2)-x(1)); 
y2 = 2.2 - 0.00008*(x(2)-x(1)); 
plot(x/1000,[y1 2.2],'k-'),hold on
plot(x/1000,[y2 2.2],'k-')
xlim([8 42]),ylim([0 2.5])