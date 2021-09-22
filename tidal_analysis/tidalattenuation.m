%% calc tidal range attenuation:
clear all,close all,clc
load('C:\GLOVER\output\myanmar\aqd\aqd_sep19_hc.mat'),hc=aqd;
load('C:\GLOVER\output\myanmar\aqd\aqd_sep19_lc.mat'),lc=aqd;
load('C:\GLOVER\output\myanmar\aqd\aqd_sep19_ag.mat'),ag=aqd;
load('C:\GLOVER\output\myanmar\longterminst\BogaleRiverInstruments.mat'),br = BogaleRiver;


figure;
plot(hc.time,hc.depth,'k--'),hold on
plot(lc.time,lc.depth,'k-.')
plot(ag.time,ag.depth,'k:')
plot(br.datenum,br.FredaDepth,'b--')
plot(br.datenum,br.MeinmahlaDepth,'-b')
xlim([datenum('9/27/2019') datenum('10/1/2019')])