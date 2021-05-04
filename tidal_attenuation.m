%% calc tidal range at all aqds:

clear all,close all,clc

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
figure;
plot(loc(hf),tr(hf),'ko'),hold on
plot(loc(lf),tr(lf),'ro')
% plot(loc(hf),avg_range(hf),'k*')
% plot(loc(lf),avg_range(lf),'r*')
legend({'high flow','low flow'})
% legend({'hf max range','lf max range','hf avg range','lf avg range'})
% ylim([0 3.6]),xlim([0 45])
ylabel('Max. observed tidal range (m)'),xlabel('Distance from mouth (km)')