clear all,close all,clc
cd C:\GLOVER\output\myanmar\aqd

F = dir('aqd*.mat');

% figure out the good stretch of data in each aqd record
F(1).idx = 1:170;
F(2).idx = 300:450;
F(3).idx = 1:250;
F(4).idx = 1000:1700;
F(5).idx = 1300:2000;
F(6).idx = 1:1100;
F(7).idx = 2150:2900;
F(8).idx = 100:225;

% look at all good deployments together and all bad deployments together
% top row is HF (Sept)
% bottom row is LF (Mar)
% Col 1 = HC, Col 2 = LC, Col 3 = Ag
% dont bother with 8 (cyclone), plot 4 and 7 on same plot
F(1).plotloc = 6;
F(2).plotloc = 4;
F(3).plotloc = 5;
F(4).plotloc = 1;
F(5).plotloc = 2;
F(6).plotloc = 3;
F(7).plotloc = 1;

% good = [1,2,5,6,8];bad = [3,4,7];

figure; 
for jj=1:6
    load(F(jj).name)
    aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.05;
    aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
    
    if jj==6
        aqd.ssc1 = aqd.ext1*0.5;
        aqd.ssc2 = aqd.ext2*0.5;
        aqd.ssc1(aqd.ssc1<2)=NaN;aqd.ssc2(aqd.ssc2<2)=NaN;
        aqd.ssc1(3500:end)=NaN;
        aqd.ssc1=movmean(aqd.ssc1,3);aqd.ssc2=movmean(aqd.ssc2,3);
    end
    ssc_mean = nanmean([aqd.ssc1,aqd.ssc2],2);

    subplot(2,3,F(jj).plotloc)
    scatter(aqd.spd_mean(F(jj).idx),aqd.depth(F(jj).idx),[],ssc_mean(F(jj).idx))
    hold on
    colorbar,caxis([0 200]),colormap(cmocean('amp'))

    if F(jj).plotloc==1 || F(jj).plotloc==4
        axis([-1.2 1.2 3 6.5]),title('High Connectivity, Mang')
    elseif F(jj).plotloc==2 || F(jj).plotloc==5
        axis([-0.5 0.5 0 3.5]),title('Low Connectivity, Mang')
    elseif F(jj).plotloc==3 || F(jj).plotloc==6
        axis([-0.5 0.5 0 3.5]),title('Low Connectivity, Agri')
    end

end
subplot(2,3,1),ylabel('High Flow','FontWeight','Bold')
subplot(2,3,4),ylabel('Low Flow','FontWeight','Bold')


figure; 
for jj=1:3
    load(F(jj).name)
    aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.05;
    aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);

    subplot(1,3,F(jj).plotloc-3)
    scatter(aqd.spd_mean(F(jj).idx),aqd.depth(F(jj).idx),[],aqd.sal(F(jj).idx))
    colorbar,caxis([8 18]),colormap(cmocean('amp'))
    
    if F(jj).plotloc==4
        axis([-1.2 1.2 3 6.5]),title('LF High Connectivity, Mang')
    elseif F(jj).plotloc==5
        axis([-0.5 0.5 0 3.5]),title('LF Low Connectivity, Mang')
    elseif F(jj).plotloc==6
        axis([-0.5 0.5 0 3.5]),title('LF Low Connectivity, Agri')
    end

end


% 
% figure; 
% for jj=1:7
%     load(F(jj).name)
%     aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.05;
%     aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
%     ssc_mean = nanmean([aqd.ssc1,aqd.ssc2],2);
%     
%     subplot(2,3,F(jj).plotloc)
%     scatter(aqd.spd_mean,ssc_mean)
%     hold on
% 
%     if F(jj).plotloc==1 || F(jj).plotloc==4
% %         axis([-1 1 0 6.5])
%         title('High Connectivity, Mang')
%     elseif F(jj).plotloc==2 || F(jj).plotloc==5
% %         axis([-0.5 0.5 0 4])
%         title('Low Connectivity, Mang')
%     elseif F(jj).plotloc==3 || F(jj).plotloc==6
% %         axis([-0.5 0.5 0 2.5])
%         title('Low Connectivity, Agri')
%     end
% 
% end
% subplot(2,3,1),ylabel('High Flow','FontWeight','Bold')
% subplot(2,3,4),ylabel('Low Flow','FontWeight','Bold')
%%
clear all,close all,clc
load('aqd_sep19_ag.mat'),hf=aqd;
load('aqd_mar18_ag.mat'),lf=aqd;

hf.spd_mean(hf.spd_mean>0.05) = hf.spd_mean(hf.spd_mean>0.05)-0.05;
hf.spd_mean(hf.slope<0)=(-1)*hf.spd_mean(hf.slope<0);
lf.spd_mean(lf.spd_mean>0.05) = lf.spd_mean(lf.spd_mean>0.05)-0.05;
lf.spd_mean(lf.slope<0)=(-1)*lf.spd_mean(lf.slope<0);

figure;
subplot(2,2,1),pcolor(hf.time,hf.z,hf.spd')
ylim([0 3]),colorbar,shading flat,caxis([0 1])
yyaxis right,plot(hf.time,hf.spd_mean,'k'),refline(0,0)
title('high flow')
subplot(2,2,3),pcolor(hf.time,hf.z,hf.dir')
ylim([0 3]),colorbar,shading flat,caxis([0 360])
subplot(2,2,2),pcolor(lf.time,lf.z,lf.spd')
ylim([0 3]),colorbar,shading flat,caxis([0 1])
yyaxis right,plot(lf.time,lf.spd_mean,'k'),refline(0,0)
title('low flow')
subplot(2,2,4),pcolor(lf.time,lf.z,lf.dir')
ylim([0 3]),colorbar,shading flat,caxis([0 360])

%%
clear all,close all,clc
cd C:\GLOVER\output\myanmar\

% sept2019 aquadopp data
load('aqd\aqd_sep19_hc.mat')
aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.05;
aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
hc=aqd;
load('aqd\aqd_sep19_lc.mat')
aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.05;
aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
lc=aqd;
load('aqd\aqd_sep19_ag.mat')
aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.05;
aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
ag=aqd;


figure;
subplot(131)
scatter(hc.spd_mean,hc.depth_elev,[],hc.ssc1)
colorbar,caxis([0 200]),colormap(cmocean('amp'))
subplot(132)
scatter(lc.spd_mean,lc.depth_elev,[],lc.ssc1)
colorbar,caxis([0 200]),colormap(cmocean('amp'))
subplot(133)
scatter(ag.spd_mean,ag.depth_elev,[],ag.ssc1)
colorbar,caxis([0 200]),colormap(cmocean('amp'))


%% CM Spring 2021 figure: compare agri to LC
clear all,close all,clc

load('aqd_sep17_lc.mat')
aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.05;
aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
lc = aqd;
load('aqd_mar18_ag.mat')
aqd.spd_mean(aqd.spd_mean>0.05) = aqd.spd_mean(aqd.spd_mean>0.05)-0.05;
aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);
ag=aqd;clear aqd

figure; 
subplot(121)
scatter(lc.spd_mean,lc.depth,[],lc.ssc1,'.')
title('Meinmahla Island')
subplot(122)
scatter(ag.spd_mean(1:450),ag.depth(1:450),[],ag.ssc1(1:450),'.')
colorbar,title('Agricultural field')

for jj=1:2
    subplot(1,2,jj)
    caxis([0 100]),colormap(cmocean('amp'))
    axis([-0.5 0.5 0 3.5])
    ylabel('Water Depth (m)')
    xlabel('Water Velocity (m/s)')
end

%%
clear all,close all,clc

load('aqd_sep17_hc.mat')
% aqd.spd_mean(aqd.spd_mean>0.5) = aqd.spd_mean(aqd.spd_mean>0.5)-0.5;
aqd.spd_mean = aqd.spd_mean-0.2;
aqd.spd_mean(aqd.slope<0)=(-1)*aqd.spd_mean(aqd.slope<0);

figure; 
scatter(aqd.spd_mean,aqd.depth,[],aqd.ssc1,'.')
colorbar
caxis([0 100]),colormap(cmocean('amp'))
axis([-1.2 1.2 3 6.5])
ylabel('Water Depth (m)')
xlabel('Water Velocity (m/s)')
