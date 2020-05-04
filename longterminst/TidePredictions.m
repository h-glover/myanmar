clear all,close all,clc
load('PatheinRiverInstruments.mat')

% Only use first deployment:
idx=find(PatheinRiver.datenum>datenum('27-Jun-2017') &...
    PatheinRiver.datenum<datenum('20-Jan-2018'));

tvec=PatheinRiver.datenum(idx)';
wl=PatheinRiver.Depth(idx);

coef=ut_solv(tvec-6.5/24,wl,[],16,'auto','NoTrend');
PatheinRiver.ut_Depth=ut_reconstr(PatheinRiver.datenum-6.5/24,coef)';
save('PatheinRiverInstruments.mat','PatheinRiver')
% figure;
% plot(PatheinRiver.datenum,PatheinRiver.Depth,'k')
% hold on
% plot(PatheinRiver.datenum,wl_pred,'r--');

%%
clear all,close all,clc
load('YangonRiverInstruments.mat')

% figure;
% plot(YangonRiver.datenum,YangonRiver.Depth)

% Only use second deployment:
idx=find(YangonRiver.datenum>datenum('17-Jul-2018') &...
    YangonRiver.datenum<datenum('15-Nov-2018'));

tvec=YangonRiver.datenum(idx);
wl=YangonRiver.Depth(idx);

coef=ut_solv(tvec-6.5/24,wl,[],16,'auto','NoTrend');
YangonRiver.ut_Depth=ut_reconstr(YangonRiver.datenum-6.5/24,coef)';
save('YangonRiverInstruments.mat','YangonRiver')
figure;
plot(YangonRiver.datenum,YangonRiver.Depth,'k')
hold on
plot(YangonRiver.datenum,YangonRiver.ut_Depth,'r--');

%%
clear all,close all,clc
load('BogaleRiverInstruments.mat')

% Use first few deployments for Freda dock sensor:
idx=find(BogaleRiver.datenum>datenum('13-Sep-2017') &...
    BogaleRiver.datenum<datenum('26-Jan-2019'));

tvec=BogaleRiver.datenum(idx)';
wl=BogaleRiver.FredaDepth(idx);

coef=ut_solv(tvec-6.5/24,wl,[],16,'auto','NoTrend');
BogaleRiver.ut_FredaDepth=ut_reconstr(BogaleRiver.datenum-6.5/24,coef)';
% 
% figure;
% plot(BogaleRiver.datenum,BogaleRiver.FredaDepth,'k')
% hold on
% plot(BogaleRiver.datenum,BogaleRiver.ut_FredaDepth,'r--');
% 

% Use first few deployments for Meinmahla dock sensor:
idx=[];tvec=[];wl=[];
idx=find(BogaleRiver.datenum>datenum('24-Feb-2017') &...
    BogaleRiver.datenum<datenum('11-Sep-2017'));

tvec=BogaleRiver.datenum(idx)';
wl=BogaleRiver.MeinmahlaDepth(idx);

coef=ut_solv(tvec-6.5/24,wl,[],16,'auto','NoTrend');
BogaleRiver.ut_MeinmahlaDepth=ut_reconstr(BogaleRiver.datenum-6.5/24,coef)';

% figure;
% plot(BogaleRiver.datenum,BogaleRiver.MeinmahlaDepth,'k')
% hold on
% plot(BogaleRiver.datenum,BogaleRiver.ut_MeinmahlaDepth,'r--');

save('BogaleRiverInstruments.mat','BogaleRiver')
