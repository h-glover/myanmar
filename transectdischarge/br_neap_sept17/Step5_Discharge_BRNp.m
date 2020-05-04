clear all
load BR_Sept17_neap.mat
for ii=1:length(adcp)
MeasuredDis(ii)=nansum(nansum(adcp(ii).alongComplete./100)).*0.5;

MeasuredTime(ii)=mean(time);
end

InterpTimes=days([min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.4)]);
InterpDis=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],[MeasuredDis,MeasuredDis(1)],InterpTimes,'linear');
InterpRes=mean(InterpDis);

%save('BR_neap_discharge_interp','InterpDis','InterpTimes','InterpRes')
%%
clear all
load BR_Sept17_Neap_SedFlux.mat
for ii=1:length(adcp)
MeasuredTime(ii)=mean(adcp(ii).time);
end
MeasuredDis=horzcat(adcp.MeasuredDis);
InterpTimes=days([min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.4)]);
InterpDis=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],[MeasuredDis,MeasuredDis(1)],InterpTimes,'linear');
InterpRes=mean(InterpDis);

MeasuredFlux=horzcat(adcp.sigmaMeasuredSedFlux);
InterpFlux=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],[MeasuredFlux,MeasuredFlux(1)],InterpTimes,'linear');
InterpResFlux=mean(InterpFlux);


save('BR_neap_discharge_interp','Interp*','Measured*')
%%


figure(1000)
clf
textsize=14;
set(gcf,'position',[1000 656 958 682])
ax1=tight_subplot(1,1,.1,.1,.1);
ax2=tight_subplot(1,1,.1,.1,.1);

set(gcf,'color','w')
axes(ax1(1))

scatter(MeasuredTime,MeasuredDis,50,'b')
hold on
plot(InterpTimes,InterpDis,':b','linewidth',2)
% load Xingu_2011_pres.mat
% [hAx,hLine1,hLine2]=plotyy(InterpTimes,InterpDis,time(31:end-4),pres(31:end-4)-mean(pres(31:end-4)))
datetick
set(gca,'fontsize',16)
xlabel('Time','fontsize',16)
ylabel('Discharge (m^3 s^{-1})','fontsize',16,'color','b') % left y-axis 

residual=mean(InterpDis);
plot(InterpTimes,ones(length(InterpTimes)).*residual,'--b','linewidth',2)

t=text(min(InterpTimes),residual+1000,strcat(num2str(round(residual)),' m^3 s^-^1'));
t.FontSize=16;
t.Color='b';
title('Bogale River | Spring Tide | Low Discharge','fontsize',16)
pos=get(gca,'outerposition');
grid on

% %%
% axes(ax2(1))
% plot(ChanP.time(3790:3868),ChanP.depth(3829:3907)-nanmean(ChanP.depth(3829:3907)),'r')
% set(ax2,'ylim',[-2,2])
% % set(gca,'fontsize',textsize)
% % xlabel('Date')
% % ylabel('Water level (m)','fontsize',textsize)
% % title('Yangon River Pressure Sensor','fontsize',textsize+4)
% % xlim([find(ChanP.datetime=='18-Sep-2017 06:00:00'),find(ChanP.datetime=='18-Sep-2017 21:00:00')])
% datetick
% % set(gca,'outerposition',pos)
% set(ax2,'yaxislocation','right')
% ylabel('Water level (m)','color','r')
% set(ax2,'xticklabels','')
% set(ax2,'xtick',[])
% set(ax2,'color','none')
% set(ax2,'fontsize',textsize)
% % set(ax2,'ytick',[0,20,40])
% set(ax2,'ycolor','r')
% set(ax1,'ycolor','b')
% set(ax2,'boxstyle','full')
