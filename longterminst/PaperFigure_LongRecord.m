% long term record figure
clear all,close all,clc
load('BogaleRiverInstruments.mat')
load('PatheinRiverInstruments.mat')
load('YangonRiverInstruments.mat')

BogaleRiver.MeinmahlaDepth=BogaleRiver.MeinmahlaDepth-nanmean(BogaleRiver.MeinmahlaDepth);
BogaleRiver.FredaDepth=BogaleRiver.FredaDepth-nanmean(BogaleRiver.FredaDepth);
PatheinRiver.Depth=PatheinRiver.Depth-nanmean(PatheinRiver.Depth);

BogaleRiver.FredaSSC=movmean(BogaleRiver.FredaSSC,35,'omitnan');
BogaleRiver.MeinmahlaSSC=movmean(BogaleRiver.MeinmahlaSSC,35,'omitnan');
PatheinRiver.SSC=movmean(PatheinRiver.SSC,30,'omitnan');
idx=find(BogaleRiver.datenum<datenum('Sept-12-2017'));
figure;
subplot(311)
plot(YangonRiver.datenum,YangonRiver.Depth,'b'),hold on
plot(BogaleRiver.datenum,BogaleRiver.FredaDepth,'k')
plot(BogaleRiver.datenum(idx),BogaleRiver.MeinmahlaDepth(idx),'k')
plot(PatheinRiver.datenum,PatheinRiver.Depth,'r')
legend('Yangon','Bogale','','Pathein')
ylabel('water level (m)'),ylim([-3.6 3.6])

subplot(312)
plot(BogaleRiver.datenum,BogaleRiver.FredaSSC,'k'),hold on
plot(BogaleRiver.datenum(62000:end),BogaleRiver.MeinmahlaSSC(62000:end),'k')
plot(PatheinRiver.datenum,PatheinRiver.SSC,'r')
ylim([0 1200])
ylabel('SSC (mg/L)')

subplot(313)
plot(BogaleRiver.datenum,BogaleRiver.MeinmahlaSal,'k'),hold on
plot(YangonRiver.datenum,YangonRiver.Sal,'b')
ylabel('Salinity')
%%
for jj=1:3
    subplot(3,1,jj)
    
    xlim([datenum('Feb-01-2017') datenum('Dec-31-2018')])
    datetick('x','mmm ''yy','keeplimits')
end

%% PLOT AYE TOTAL DISCHARGE



clear all,close all,clc

fid=fopen('Furuichi_T2.csv');
T=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter',',');
fclose (fid);

disch=T{2};
for jj=2:12
	disch=horzcat(disch,T{jj});
end
disch=[disch(:,end),disch,disch(:,1)];
disch_std=std(disch(1:31,:));
upp=disch(32,:)+disch_std;
dwn=disch(32,:)-disch_std;
x=datenum('15-Dec-2018'):datenum(0,0,30,0,0,0):datenum('15-Jan-2020');

figure;
plot(x,disch(32,:),'k'),hold on
plot(x,upp,'r')
plot(x,dwn,'r')
datetick('x','mmm','keeplimits')
xlim([datenum('25-Dec-2018') datenum('25-Dec-2019')])