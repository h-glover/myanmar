clear all%,close all,clc
%load('YR_Sept17_SedFlux.mat')
load('YR_Sept17_discharge_interp.mat')

%
F=figure;
F.Position=[10 50 1000 500];
yyaxis right
plot(MeasuredTime,MeasuredFlux,'r*')
hold on
plot(InterpTimes,InterpFlux,'r--')
r=refline(0,InterpResFlux); r.Color='r';r.LineStyle=':';
ylim([-7 7])

%
ylabel('Sed flux, t/s')
yyaxis left
plot(InterpTimes,InterpDis,'k--')
plot(MeasuredTime,MeasuredDis,'k-',MeasuredTime,MeasuredDis,'ko')
r=refline(0,InterpRes); r.Color='k';r.LineStyle=':';
ax=gca;ax.YColor='k';
%sigmaMeasuredDis=vertcat(adcp.sigmaMeasuredDis);
%plot(MeasuredTime,sigmaMeasuredDis,'b*')
ylim([-18000 18000]),ylabel('Discharge, m^3/s')
xlim([InterpTimes(1)-1/48 InterpTimes(end)+1/48])
datetick('x','HH:MM','keeplimits')
xlabel(datestr(floor(MeasuredTime(1))))
title('Yangon River: High Flow: Spring')
t=text(min(InterpTimes),InterpRes+1000,strcat(num2str(round(InterpRes)),' m^3 s^-^1'));
t.Color='k';
yyaxis right
t=text(MeasuredTime(end-5),InterpResFlux-1,strcat(num2str(round(InterpResFlux,2)),' t/s'));
t.Color='r';
