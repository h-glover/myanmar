% Process aquadopp data from MeinMaHla Island
clear all,close all,clc
aqd=aqd_process('BigCh01',0.05,2,2,[],[],0.5);

% remove air pressure signal using Met station and add HAB:
load('BogaleRiverWeather.mat')
Weather.AtmPres=fillmissing(Weather.AtmPres,'nearest');
airpres=interp1(Weather.datenum,Weather.AtmPres,aqd.time)./100; % in dbar
figure;
plot(aqd.pres,'k'),hold on,
plot(airpres,'r')
nanmean(aqd.pres(([1:53,601:end]))-airpres(([1:53,601:end])))
