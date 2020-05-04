% aqd plotting
clear all%,close all,clc
load('MMsm_Mar18_aqd')
smch=aqd;
load('MMLG_Mar18_aqd')
lgch=aqd;
load('MMAG_Mar18_aqd')
agch=aqd; clear aqd

smch.temp(smch.pres<0.2)=NaN;
smch.ext1(smch.pres<0.2)=NaN;
smch.ext2(smch.pres<0.2)=NaN;
smch.sal(smch.pres<0.2 | smch.sal<10)=NaN;
smch.pres(smch.pres<0.2)=NaN;

lgch.temp(lgch.pres<2)=NaN;
lgch.ext1(lgch.pres<2)=NaN;
lgch.ext2(lgch.pres<2)=NaN;
lgch.sal(lgch.pres<2)=NaN;
lgch.pres(lgch.pres<0.2)=NaN;

agch.temp(agch.pres<0.2)=NaN;
agch.ext1(agch.pres<0.2 | agch.ext1<140)=NaN;
agch.ext2(agch.pres<0.2)=NaN;
agch.sal(agch.pres<0.2)=NaN;
agch.pres(agch.pres<0.2)=NaN;


lgch.ssc1=lgch.ext1.*0.0688 + 50;
smch.ssc1=smch.ext1.*0.0688 + 50;
agch.ssc1=agch.ext1.*0.0688 + 50;


lgch.slope=lgch.pres(2:end)-lgch.pres(1:end-1);
lgch.slope=[lgch.slope(1);lgch.slope]*12;
smch.slope=smch.pres(2:end)-smch.pres(1:end-1);
smch.slope=[smch.slope(1);smch.slope]*12;
agch.slope=agch.pres(2:end)-agch.pres(1:end-1);
agch.slope=[agch.slope(1);agch.slope]*12;


%% plot turb/sal/pres as a time series

figure;
subplot(411)
plot(smch.time,smch.pres,'b'),hold on
plot(lgch.time,lgch.pres,'k')
%plot(agch.time,agch.pres,'Color',[0,0.7,0])
ylim([0 7]),ylabel('Depth (m)')
xlim([agch.time(2) smch.time(end)])
datetick('x','dd','keeplimits')
legend('Dead-end Channel','Main Channel','Agri Channel')

subplot(412)
plot(smch.time,smch.temp,'b'),hold on
plot(lgch.time,lgch.temp,'k')
%plot(agch.time,agch.temp,'Color',[0,0.7,0])
xlim([agch.time(2) smch.time(end)])
ylim([27 32]),ylabel('Temp (C)')
datetick('x','dd','keeplimits')

subplot(413)
plot(smch.time,smch.sal,'b'),hold on
plot(lgch.time,lgch.sal,'k')
%plot(agch.time,agch.sal,'Color',[0,0.7,0])
xlim([agch.time(2) smch.time(end)])
ylim([0 20]),ylabel('Salinity')
datetick('x','dd','keeplimits')

subplot(414)
plot(smch.time,smch.ssc1,'b'),hold on
plot(lgch.time,lgch.ssc1,'k')
%plot(agch.time,agch.ssc1,'Color',[0,0.7,0])
xlim([agch.time(2) smch.time(end)])
datetick('x','dd','keeplimits')
ylim([0 4000]),ylabel('SSC (mg/L) at 20 cmab')
xlabel(['Day in ',datestr(lgch.time(1),'mmm ''yy')])

%% slope-salinity plots colored by ssc
figure;
subplot(2,2,3)
scatter(lgch.slope,lgch.sal,[],lgch.ssc1,'.')
caxis([0 700])
C=[0.8*ones([length(lgch.time),1]),linspace(0.8,0,length(lgch.time))',...
    linspace(0.8,0,length(lgch.time))'];
colormap(C)
colorbar
axis([-1 1 12 19])
title('Low Flow: Large Channel')
xlabel('Slope (m/hr)'),ylabel('Salinity')
subplot(224)
scatter(smch.slope,smch.sal,[],smch.ssc1,'.')
colormap(C)
colorbar,caxis([0 700])
axis([-1 1 12 19])
title('Low Flow: Dead-end Channel')
xlabel('Slope (m/hr)')

%% salinity ssc plots colored by slope
C1=flipud([ones([30,1]),linspace(1,0,30)',linspace(1,0,30)']);
C2=[linspace(1,0,30)',linspace(1,0,30)',ones([30,1])];
C=[C1(1:end-4,:);C2(5:end,:)];


subplot(2,2,3)
scatter(lgch.sal,lgch.ssc1,[],lgch.slope,'.')
caxis([-1 1]),colorbar,colormap(C)
%axis([12 19 0 2000])
title('High Flow: Large Channel')
ylabel('SSC')
subplot(2,2,4)
scatter(smch.sal,smch.ssc1,[],smch.slope,'.')
caxis([-1 1]),colorbar
%axis([12 19 0 2000])
title('High Flow: Dead-end Channel')

%% all three locations on slope, sal, ssc plots
figure;
subplot(131)
scatter(lgch.slope,lgch.sal,[],lgch.ssc1,'.')
C=[0.8*ones([length(lgch.time),1]),linspace(0.8,0,length(lgch.time))',...
    linspace(0.8,0,length(lgch.time))'];
colormap(C)
colorbar%caxis([0 700])
axis([-1 1 12 19])
title('Low Flow: Large Channel')
xlabel('Slope (m/hr)'),ylabel('Salinity')
%
subplot(132)
scatter(smch.slope,smch.sal,[],smch.ssc1,'.')
colormap(C)
colorbar,%caxis([0 700])
axis([-1 1 12 19])
title('Low Flow: Dead-end Channel')
xlabel('Slope (m/hr)')
%
subplot(133)
figure;
scatter(agch.slope,agch.sal,[],agch.ssc1,'.')
colormap(C)
colorbar,caxis([0 800])
axis([-1 1 2 10])
title('Low Flow: Agri Channel')
xlabel('Slope (m/hr)')


%% all three locations on slope, sal, ssc plots
C=[0.8*ones([length(lgch.time),1]),linspace(0.8,0,length(lgch.time))',...
    linspace(0.8,0,length(lgch.time))'];

figure;
scatter(lgch.slope,lgch.sal,[],lgch.ssc1,'x')
hold on
scatter(smch.slope,smch.sal,[],smch.ssc1,'o')
scatter(agch.slope,agch.sal,[],agch.ssc1,'d')
colormap(C)
colorbar,caxis([0 600])
axis([-1 1 2 18.5])
legend('HC','LC','AG')
xlabel('Slope (m/hr)'),ylabel('Salinity')