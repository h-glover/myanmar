clear all,close all,clc
%load('BR_Sept17_Neap_SedFlux.mat')
load('BR_neap_discharge_interp.mat')

%
F=figure;
F.Position=[10 50 1000 500];
yyaxis right
plot(MeasuredTime,MeasuredFlux,'r*')
hold on
plot(InterpTimes,InterpFlux,'r--')
r=refline(0,InterpResFlux); r.Color='r';r.LineStyle=':';
ylim([-7 7])
hold on
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
title('Bogale River: High Flow: Neap')
t=text(max(MeasuredTime),InterpRes+1000,strcat(num2str(round(InterpRes)),' m^3 s^-^1'));
t.Color='k';
yyaxis right
t=text(max(MeasuredTime),InterpResFlux-1,strcat(num2str(round(InterpResFlux,2)),' t/s'));
t.Color='r';
