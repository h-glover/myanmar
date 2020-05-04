% clear all,clc
% 
% aqd=aqd_process('BigCh01',0.05,2,2,[],[],0.5);
% 
% aqd.depth=NaN(length(aqd.time),1);
% airoffset=nanmean(aqd.pres([1:60,end-10:end]));
% aqd.depth=aqd.pres-airoffset;
% aqd.depth(aqd.depth<1)=NaN;
% 
% save('MainCh_Sep19','aqd')
% clear all,clc
% 
% aqd=aqd_process('BigCh01',0.05,2,2,[],[],0.5);
% 
% aqd.depth=NaN(length(aqd.time),1);
% airoffset=nanmean(aqd.pres([1:60,end-10:end]));
% aqd.depth=aqd.pres-airoffset;
% aqd.depth(aqd.depth<1)=NaN;
% 
% save('MainCh_Sep19','aqd')

%% survey map data
% clear all,close all,clc
% 
% fid=fopen('SurveyPts2019.csv');
% T=textscan(fid,'%*f %*f %*f %s %f %f %f %f %f',...
%     'Delimiter',',','HeaderLines',1);
% fclose(fid);
% 
% 
% survey.site=T{1};
% survey.Lat=T{3};
% survey.Lon=T{2};
% survey.Z=T{4};
% survey.Zcorr=T{5};
% survey.Zref=T{6};
% save('SurveyPts2019','survey')

%% process data
clear all,close all,clc

% compare elevation to water level
load('SurveyPts2019.mat')
load('C:\GLOVER\Myanmar\AQD\BR_Aqd_Sept19\MainCh_Sep19.mat')
load('C:\GLOVER\Myanmar\LongTermInst\BogaleRiverInstruments.mat')
% figure;
% subplot(311),plot(aqd.head)
% subplot(312),plot(aqd.depth)
% subplot(313),plot(aqd.pitch,'k'),hold on,plot(aqd.roll,'r')

tvec=datenum('28-Sept-19 08:00'):datenum(0,0,0,0,1,0):datenum('28-Sept-19 09:00');
wl=interp1(aqd.time,aqd.depth,tvec);
wl(tvec==datenum('28-Sept-19 08:34'))

% 0834 water level, BR36=Row166, depth=4.5686, Zraw=-43.1357;
% next find mean water level compared to this...

figure;
plot(BogaleRiver.datenum,BogaleRiver.FredaDepth,'k'),hold on
plot(aqd.time,aqd.depth,'b')
plot(tvec,wl,'r')
xlim([datenum('20-Sep-19') datenum('30-Sep-19')])

