clear all%,close all,clc
%load('BR_Mar18_neap_SedFlux.mat')
load('BR_Mar18_Neap_discharge_interp')
%load('BR_neap_Salinity')
%
F=figure;
F.Position=[10 50 1000 500];
yyaxis right
plot(MeasuredTime,MeasuredFlux,'r*')
hold on
plot(InterpTimes,InterpFlux,'r--')
ylim([-4 4])
r=refline(0,InterpResFlux); r.Color='r';r.LineStyle=':';
%
ylabel('Sed flux, t/s')
yyaxis left
plot(InterpTimes,InterpDis,'k--')
plot(MeasuredTime,MeasuredDis,'k-',MeasuredTime,MeasuredDis,'ko')
r=refline(0,InterpRes); r.Color='k';r.LineStyle=':';
ax=gca; ax.YColor='k';
%sigmaMeasuredDis=vertcat(adcp.sigmaMeasuredDis);
%plot(MeasuredTime,sigmaMeasuredDis,'b*')
ylim([-11000 11000]),ylabel('Discharge, m^3/s')
xlim([InterpTimes(1)-1/48 InterpTimes(end)+1/48])
datetick('x','HH:MM','keeplimits')
xlabel(datestr(floor(MeasuredTime(1))))
title('Bogale River: Low Flow: Neap')
t=text(min(InterpTimes),InterpRes+1000,strcat(num2str(round(InterpRes)),' m^3 s^-^1'));
t.Color='k';
yyaxis right
t=text(max(MeasuredTime),InterpResFlux-1,strcat(num2str(round(InterpResFlux,2)),' t/s'));
t.Color='r';

%%
idx=find(stnidx==2);
ax2 = axes('Position',ax.Position,'Color','none','Box','off');
ax2.Position(4)=0.06;ax2.Position(2)=0.9;
p=pcolor(ax2,timevec(idx),depthvec,salvec(:,idx));
shading interp,axis ij
p.EdgeColor='interp';
c=colorbar;caxis([0 12])
c.Label.String='Sal';
c.Location='eastoutside';c.Position(4)=0.08;c.Position(2)=0.89;
ylim([2 4])
xlim([InterpTimes(1)-1/48 InterpTimes(end)+1/48])
datetick('x','HH:MM','keeplimits')
ax2.YTickLabels=[];ax2.XTickLabels=[];

