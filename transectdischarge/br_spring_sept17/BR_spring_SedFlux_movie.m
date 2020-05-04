% dsfs
clear all,close all,clc

load('BR_Mar18_Spring_SedFlux.mat')
crosssection=linspace(0,1,1000);

v=VideoWriter('BR_Spring_LowFlow_TransectsIS');
v.FrameRate=1;
open(v)
for jj=1:length(adcp)
    fig=figure('visible','off');
    subplot(311)
    pcolor(0:1000,crosssection,adcp(jj).sigmaAlongComplete.*100),shading flat,axis ij
    c=colorbar;caxis([-2 2]),ylabel('Discharge')
    c.Label.String='m^3/s';
    hold on
    plot(adcp(jj).CTDlocations,zeros([1,length(adcp(jj).CTDlocations)]),'rv')
    title(['Transect ',num2str(jj),' ',datestr(adcp(jj).time(1),'mm/dd/yy HH:MM')])
    subplot(312)
    pcolor(0:1000,crosssection,adcp(jj).ssc),shading flat,axis ij
    c=colorbar;caxis([0 1200]),ylabel('SSC')
    c.Label.String='mg/L';
    hold on
    plot(adcp(jj).CTDlocations,zeros([1,length(adcp(jj).CTDlocations)]),'rv')
    subplot(313)
    pcolor(0:1000,crosssection,adcp(jj).sigmaSedFlux./1000),shading flat,axis ij
    c=colorbar;caxis([-12 12]),ylabel('Sed Flux')
    c.Label.String='g/s';
    hold on
    plot(adcp(jj).CTDlocations,zeros([1,length(adcp(jj).CTDlocations)]),'rv')
    F(jj)=getframe(fig,[0 0 560 420]);
    writeVideo(v,F(jj))
    
end
close(v)

%% Look at CTD profiles individually
clear all,close all,clc

load('BR_Mar18_Spring_SSC.mat')
alldepth=horzcat(ssc.Depth);
allssc=horzcat(ssc.SSC);
idx=1:3:length(ssc)*3;
for jj=idx
    figure;
    plot(allssc(:,jj),alldepth(:,jj),'k')
    hold on
    plot(allssc(:,jj+1),alldepth(:,jj+1),'r')
    plot(allssc(:,jj+1),alldepth(:,jj+1),'b')
    axis ij,xlim([0 1500]),ylim([0 14])
end



