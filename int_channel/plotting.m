%
clear all,close all,clc

cd C:\GLOVER\output\myanmar\aqd

F = dir('*.mat');
figure;
for jj=1:length(F)
    load(F(jj).name)
    
    aqd.spd_mean(aqd.slope<0) = (-1)*aqd.spd_mean(aqd.slope<0);
    vel = aqd.spd(:,2)-0.05;
    vel(aqd.slope<0) = (-1)*vel(aqd.slope<0);
    
    figure(1)
    subplot(2,4,jj)
    scatter(vel,aqd.depth,[],aqd.ssc1,'.'),hold on
    scatter(vel,aqd.depth,[],aqd.ssc2,'.')
    xlabel('velocity'),ylabel('water depth')
    caxis([0 100]),colormap('jet')
    xlim([-2 2]),ylim([0 6.5])
    title(F(jj).name(5:end-4),'Interpreter','none')
    
    figure(2)
    subplot(2,4,jj)
    scatter(aqd.slope,aqd.depth,[],aqd.ssc1,'.'),hold on
    scatter(aqd.slope,aqd.depth,[],aqd.ssc2,'.')
    caxis([0 100]),colormap('jet')
    xlim([-1 1]),ylim([0 6.5])
    xlabel('slope (m/hr)'),ylabel('water depth')
    title(F(jj).name(5:end-4),'Interpreter','none')

end
figure(1),colorbar
figure(2),colorbar

%% scatter plots with slope and ssc:

clear all,close all,clc

cd C:\GLOVER\output\myanmar\aqd

F = dir('*.mat');
for jj=1:length(F)
    loc = F(jj).name(end-5:end-4);
    
    if strcmp(loc,'ag')==1
        clr{jj} = 'g';
    elseif strcmp(loc,'lc')==1
        clr{jj} = 'b';
    elseif strcmp(loc,'hc')==1
        clr{jj} = 'k';
    end
end

figure;
subplot(121)
for jj=1:3
    load(F(jj).name)
    scatter(aqd.slope,aqd.ssc1,[],'.','MarkerEdgeColor',clr{jj}),hold on
end
title('Low flow'),xlabel('slope (m/hr)'),ylabel('ssc (mg/l)')
xlim([-1.2 1.2]),ylim([0 400])
legend({'agri','high conn','low conn'})
subplot(122)
for jj=4:6
    load(F(jj).name)
    scatter(aqd.slope,aqd.ssc1,[],'.','MarkerEdgeColor',clr{jj}),hold on
end
title('High flow'),xlabel('slope (m/hr)')
xlim([-1.2 1.2]),ylim([0 400])


