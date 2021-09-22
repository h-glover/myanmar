clear all, close all, clc

cd C:\GLOVER\output\myanmar\longterminst


load('BogaleRiverInstruments')
load('BogaleRiverWeather')
load('BogaleRiverWell')

%% river data
ylab = {'Depth (m)','Temp (C)','Salinity','SSC (mg/L)'};
figure;
subplot(411)
plot(BogaleRiver.datenum,BogaleRiver.FredaDepth,'k'),hold on
plot(BogaleRiver.datenum,BogaleRiver.MeinmahlaDepth,'r')
subplot(412)
plot(BogaleRiver.datenum,BogaleRiver.FredaTemp,'k'),hold on
plot(BogaleRiver.datenum,BogaleRiver.MeinmahlaSalTemp,'r')
subplot(413)
plot(BogaleRiver.datenum,BogaleRiver.MeinmahlaSal,'r')
subplot(414)
plot(BogaleRiver.datenum,BogaleRiver.FredaSSC,'k'),hold on
plot(BogaleRiver.datenum,BogaleRiver.MeinmahlaSSC,'r')
for jj=1:4
    subplot(4,1,jj)
    datetick('x','mm ''yy','keeplimits')
    ylabel(ylab{jj})
end

%% well data
figure;
subplot(211)
plot(BogaleWell.datenum,BogaleWell.MeinmahlaWellDepth,'k'),hold on
subplot(212)
plot(BogaleWell.datenum,BogaleWell.MeinmahlaWellTemp,'k'),hold on
for jj=1:2
    subplot(2,1,jj)
    datetick('x','mm ''yy','keeplimits')
    ylabel(ylab{jj})
end

%% weather

figure;
subplot(511)
scatter(Weather.datenum,Weather.Wind,[],Weather.WindDir,'.')
ax=gca; ax.Colormap=cmocean('phase');colorbar
ylabel('Wind Speed and Dir')
subplot(512)
plot(Weather.datenum,Weather.AtmPres)
ylabel('Atm Pres')
subplot(513)
plot(Weather.datenum,Weather.AirTemp)
ylabel('Temp (C)')
subplot(514)
yyaxis left,plot(Weather.datenum,Weather.Hum)
ylabel('Humidity (%)')
yyaxis right,plot(Weather.datenum,Weather.Rain)
ylabel('Rain (mm)')
subplot(515)
plot(Weather.datenum,Weather.PAR)
ylabel('PAR')

for jj=1:5
    subplot(5,1,jj)
    datetick('x','mm ''yy','keeplimits')
%     ylabel(ylab{jj})
end
