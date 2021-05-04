% convert pressure to depth, convert conductivity to salinity, convert
% turbidity to ssc and clean data for all time series

clear all,close all,clc

% Bogale first:
load('BogaleRiverInstruments.mat')
load('BogaleRiverWeather')

% remove atm pressure and convert to meters, AtmPres=mbar,
% MeinmahlaPres=kpa
BogaleRiver.MeinmahlaDepth=...
    (BogaleRiver.MeinmahlaPres-(Weather.AtmPres./10)).*0.101972;
BogaleRiver.FredaDepth=...
    (BogaleRiver.FredaPres-(Weather.AtmPres./10)).*0.101972;
% (future step is align the records then subtract mean depth to
% get "tidal range")

% calculate conductivity
pres_dbar=(BogaleRiver.MeinmahlaPres-(Weather.AtmPres./10))./10;
pres_dbar=nanmean(pres_dbar); 
pres_dbar=repmat(pres_dbar,1,length(BogaleRiver.MeinmahlaCond));
temp_C=BogaleRiver.MeinmahlaSalTemp;
temp_C(temp_C>35)=35;
BogaleRiver.MeinmahlaSal=conduc2sali(BogaleRiver.MeinmahlaCond./1000,temp_C,...
    pres_dbar);

% convert turb to ssc
BogaleRiver.MeinmahlaSSC=BogaleRiver.MeinmahlaTurb.*0.726;
BogaleRiver.FredaSSC=BogaleRiver.FredaTurb.*0.726;

save('BogaleRiverInstruments','BogaleRiver')
%% Clean Bogale data (only clean Depth, Salinity, SSC - not raw vals)
clear all,close all,clc

% Bogale first:
load('BogaleRiverInstruments.mat')

t=BogaleRiver;

% fix Salinity/SSC by removing low water, outliers, and smoothing
t.FredaSSC(t.FredaSSC<30)=NaN;
t.MeinmahlaSSC(t.MeinmahlaSSC<30)=NaN;
t.MeinmahlaSal(t.MeinmahlaSal<0.5)=NaN;
t.FredaSSC(t.FredaDepth<0.2)=NaN;
t.MeinmahlaSSC(t.MeinmahlaDepth<0.2)=NaN;
t.MeinmahlaSal(t.MeinmahlaDepth<0.2)=NaN;
t.MeinmahlaSal([14523:16279,53520:end])=NaN; %bad vals...mud?

nn=7;
c=repmat(1/nn,[1,nn]);
t.FredaSSC=conv(t.FredaSSC,c,'same');
t.MeinmahlaSSC=conv(t.MeinmahlaSSC,c,'same');
t.MeinmahlaSal=conv(t.MeinmahlaSal,c,'same');
figure; plot(t.MeinmahlaSSC)
figure; plot(t.FredaSSC)
figure; plot(t.MeinmahlaSal)

% fix Freda water level record by removing vals<5cm
t.FredaDepth(t.FredaDepth<0.05)=NaN;

% fix MMI water level by subtracting the difference between the mean depth
% for each deployument and the orig. mean depth (use ttide to fill missing
% low water values)
idx=[1,21068,32354,36491,65556,72099,75330,87089,132726,length(t.datenum)];
lat=16;
coef=ut_solv(t.datenum(8000:20000)',t.MeinmahlaDepth(8000:20000),[],lat,'auto');
pred=ut_reconstr(t.datenum(8000:20000),coef);
% plot(t.datenum(8000:20000),pred,'k')

runmn(1)=nanmean(pred);
baseline(1,1:idx(2))=runmn;
for jj=2:length(idx)-1
    runmn(jj)=nanmean(t.MeinmahlaDepth(idx(jj):idx(jj+1)));
    baseline=[baseline,repmat(runmn(jj),[1,(idx(jj+1)-idx(jj))])];
end
t.MeinmahlaDepth(t.MeinmahlaDepth<0.01)=NaN;
fixval=t.MeinmahlaDepth-(baseline-baseline(1));
% plot(t.datenum,fixval,'b')

% % adjust both depths to be tidal range (mean is zero)
% t.MeinmahlaDepth=t.MeinmahlaDepth-nanmean(t.MeinmahlaDepth);
% t.FredaDepth=t.FredaDepth-nanmean(t.FredaDepth);

figure;plot(t.datenum,t.FredaDepth)
figure;plot(t.datenum,t.MeinmahlaDepth)
BogaleRiver=t;
save('BogaleRiverInstruments.mat','BogaleRiver')


%% Next Pathein:

clear all,close all,clc

load('PatheinRiverInstruments.mat')
load('BogaleRiverWeather')

% remove atm pressure and convert to meters
PatheinRiver.Depth=...
    (PatheinRiver.Pres-(Weather.AtmPres./10)).*0.101972;
PatheinRiver.Depth(PatheinRiver.Depth<1)=NaN;
figure;plot(PatheinRiver.datenum,PatheinRiver.Depth),datetick

PatheinRiver.SSC=PatheinRiver.Turb.*0.726;
PatheinRiver.SSC(PatheinRiver.SSC<30)=NaN;
save('PatheinRiverInstruments','PatheinRiver')


%% Fix Yangon Depth and Salinity:
clear all,close all,clc

load('YangonRiverInstruments.mat')
%YangonRiver.Depth(38954:42841)=YangonRiver.Depth(38954:42841)+5.3;
L=length(YangonRiver.Depth)-1;
YangonRiver.Depth(1:L/2)=YangonRiver.Depth(1:L/2)-nanmean(YangonRiver.Depth(1:L/2));
YangonRiver.Depth(L/2:end)=YangonRiver.Depth(L/2:end)-nanmean(YangonRiver.Depth(L/2:end));


YangonRiver.Depth(YangonRiver.Depth>4.8  | YangonRiver.Depth<-3)=NaN;
plot(YangonRiver.Depth)
% hold on, plot(YangonRiver.Depth(1:L/2),'b')
save('YangonRiverInstruments','YangonRiver')


%% add pdf data record to MM turbidity

clear all,close all,clc
fid=fopen('test.csv');
T=textscan(fid,'%f %f','HeaderLines',1,'Delimiter',',');fclose(fid);

[xx,ia,~]=unique(T{1});
yy=T{2}(ia);

xinterp=0:0.1:648;
yinterp=interp1(xx,yy,xinterp);
xdate=linspace(datenum('8-Mar-2018'),datenum('17-Aug-2018'),length(xinterp));
xdate_interp=datenum('8-Mar-2018'):datenum(0,0,0,0,1,0):datenum('17-Aug-2018');
yinterp=interp1(xdate,yinterp,xdate_interp);

load('BogaleRiverInstruments.mat')
[~,iold,inew]=intersect(BogaleRiver.datenum,xdate_interp);
BogaleRiver.MeinmahlaTurb(iold)=yinterp(inew);

figure;
plot(BogaleRiver.datenum,BogaleRiver.MeinmahlaTurb,'k'),hold on
plot(xdate_interp,yinterp,'r--')

BogaleRiver.MeinmahlaSSC(iold)=BogaleRiver.MeinmahlaTurb(iold).*0.726;
figure;
plot(BogaleRiver.MeinmahlaSSC)
save('BogaleRiverInstruments.mat','BogaleRiver')


%% calc the water level referenced to mean water level during HF 2019
clear all,close all,clc

cd C:\GLOVER\output\myanmar\longterminst
load('BogaleRiverInstruments'),br=BogaleRiver;

msl_hf_m = nanmean(br.MeinmahlaDepth(132900:144100));
br.MeinmahlaWaterLevel=br.MeinmahlaDepth-msl_hf_m;
hhw_hf = max(br.MeinmahlaWaterLevel(132900:144100));
llw_hf = min(br.MeinmahlaWaterLevel(132900:144100));

idx = 133387:143954;
coefs = ut_solv(br.datenum(idx)',br.FredaDepth(idx),[],16,'auto');
pred = ut_reconstr(br.datenum(idx),coefs);
msl_hf_freda = nanmean(pred);
br.FredaWaterLevel = br.FredaDepth-msl_hf_freda;


BogaleRiver=br;
save('BogaleRiverInstruments','BogaleRiver')


%% fix well atm record
clear all,close all,clc

cd C:\GLOVER\output\myanmar\longterminst

load('BogaleRiverWell'),bw=BogaleWell;
load('BogaleRiverWeather')

[~,~,ib]=intersect(bw.datenum,Weather.datenum);
bw.MeinmahlaWellPres = bw.MeinmahlaWellPres - Weather.AtmPres(ib)/10;
bw.MeinmahlaWellPres(bw.MeinmahlaWellPres<0)=NaN;
bw.MeinmahlaWellDepth = bw.MeinmahlaWellPres*0.101972;
% 
% figure;
% plot(bw.MeinmahlaWellDepth)
BogaleWell = bw;

save('BogaleRiverWell.mat','BogaleWell')