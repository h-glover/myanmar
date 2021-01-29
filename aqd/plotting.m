%
clear all,close all,clc

cd C:\GLOVER\output\myanmar\aqd

F = dir('*.mat');

for jj=1:length(F)
    load(F(jj).name)
    
    
    aqd.spd_mean(aqd.slope<0) = (-1)*aqd.spd_mean(aqd.slope<0);
    vel = aqd.spd(:,2)-0.05;
    vel(aqd.slope<0) = (-1)*vel(aqd.slope<0);
    
    figure;
    subplot(121)
    scatter(vel,aqd.depth,[],aqd.ssc1),hold on
    scatter(vel,aqd.depth,[],aqd.ssc2,'.')
    xlabel('velocity'),ylabel('water depth')
    caxis([0 100])
    xlim([-2.5 2.5]),ylim([0 6.5])
    legend({'ssc2','ssc1'}),title(F(jj).name)
    subplot(122)
    scatter(aqd.slope,aqd.depth,[],aqd.ssc1),hold on
    scatter(aqd.slope,aqd.depth,[],aqd.ssc2,'.')
    colorbar,caxis([0 100])
    xlim([-1 1]),ylim([0 6.5])
    xlabel('slope (m/hr)')
end