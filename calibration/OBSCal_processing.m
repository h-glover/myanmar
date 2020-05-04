clear all,close all,clc

stime=datenum({'22-Mar-2019 13:30:02';'22-Mar-2019 13:33:30';...
    '22-Mar-2019 13:37:16';'22-Mar-2019 13:39:52';'22-Mar-2019 13:42:03';...
    '22-Mar-2019 13:44:34';'22-Mar-2019 13:46:48';'22-Mar-2019 13:51:22'});

% calculate the voltage for ext 1 and ext2 from AQD11782
fname='OBSCal01.sen';
% LOADAQDSEN Private function to return data from Aquadopp *.sen file
% [TIME VOLTAGE SOUNDSPEED HEADING PITCH ROLL PRESURE TEMP] =
% LOADAQDSEN(FNAME)
% 20091120 nowacki added external voltage readin
sen = load(fname);

% Timestamp - columns 1-6 of SEN as [month day year hour min sec]
% datenum uses [year month day hour min sec]
t = datenum(sen(:,[3 1 2 4:6]));

ext1  = sen(:,16); % external voltage input, counts (0-65535) added nowacki
ext2  = sen(:,17); % external voltage input, counts (0-65535) added nowacki

for jj=1:length(stime)
    idx=find(t>stime(jj) & t<(stime(jj)+ 25/86400));
    ext1s(jj,1)=nanmean(ext1(idx));
    ext2s(jj,1)=nanmean(ext2(idx));
end


figure;
subplot(211)
plot(t,ext1),hold on
plot(stime,ext1s,'*'),ylim([0 15000])
datetick('x','MM:SS','keeplimits')
subplot(212)
plot(t,ext2),hold on
plot(stime,ext2s,'*'),ylim([0 15000])
datetick('x','HH:MM')