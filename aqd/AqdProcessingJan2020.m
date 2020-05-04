% Process aquadopp data from MeinMaHla Island
clear all,close all,clc
aqd=aqd_process('BRmain01',0.05,2,2,[],[],0.5);
% offset for 11782 = 
offset=-9.35;

% remove bad data at beginning of deployment:
flds=fieldnames(aqd);
flds([11,19:21])=[];
for jj=1:length(flds)
    if jj<11 || jj>18
    aqd.(flds{jj})=aqd.(flds{jj})(477:end-2);
    else
    aqd.(flds{jj})=aqd.(flds{jj})(476:end-2,:);
    end
end

% load data from other platform instruments:
load('065625_20170913_1535.mat')
RBR.time=datenum(RBR.sampletimes,'yyyy-mm-dd HH:MM:SS.FFF');
aqd.rbr.depth=interp1(RBR.time,RBR.data(:,3),aqd.time);
aqd.rbr.salinity=interp1(RBR.time,RBR.data(:,7),aqd.time);
aqd.rbr.turb=interp1(RBR.time,RBR.data(:,4),aqd.time);
aqd.rbr.hsig=interp1(RBR.time,RBR.data(:,9),aqd.time);

% remove air pressure signal using Met station and add HAB:
load('BogaleRiverWeather.mat')
Weather.AtmPres=fillmissing(Weather.AtmPres,'nearest');
airpres=interp1(Weather.datenum,Weather.AtmPres,aqd.time)./100;
aqd.depth=((aqd.pres-(airpres+offset)).*1.02)+0.10;
aqd.depth(aqd.depth<3.5)=NaN;
aqd.rbr.depth=(aqd.rbr.depth-(airpres-0.03)).*1.02 - 0.45; %why 0.45 off??
aqd.rbr.depth(aqd.rbr.depth<3.5)=NaN;

% convert turbidity to SSC
aqd.SSC1=aqd.ext1.*0.122-34; % see obs calibration pdf
aqd.SSC2=aqd.ext2.*0.122-34;
aqd.rbr.ssc=aqd.rbr.turb.*2; % approx based on other cals

figure;
plot(aqd.SSC2,'k'),hold on
plot(aqd.rbr.ssc,'r')

% save('AQD_Sep17_LG','aqd')
%%  Process aquadopp data from small channel in MMI
clear all,close all,clc

aqd=aqd_process('BRdead01',0.05,2,2,[],[],0.5);
% offset for 6106 = 
offset=-9.91;

% remove out of water data:
flds=fieldnames(aqd);
flds([11,19:21])=[];
for jj=1:length(flds)
    if jj<11 || jj>18
    aqd.(flds{jj})=aqd.(flds{jj})(159:4585);
    else
    aqd.(flds{jj})=aqd.(flds{jj})(159:4585,:);
    end
end

% convert pres to depth by removing air pres and adding HAB:
load('BogaleRiverWeather.mat')
Weather.AtmPres=fillmissing(Weather.AtmPres,'nearest');
airpres=interp1(Weather.datenum,Weather.AtmPres,aqd.time)./100;
aqd.depth=((aqd.pres-(airpres+offset)).*1.02)+0.10;

% convert turbidity to SSC
aqd.SSC1=aqd.ext1.*0.111-24; % see obs calibration pdf
aqd.SSC2=aqd.ext2.*0.066-28;
aqd.SSC1(aqd.SSC1<100 | aqd.SSC1>650)=NaN;
aqd.SSC2(aqd.SSC2<100 | aqd.SSC2>650)=NaN;
% figure;
% plot(aqd.depth)
% refline(0,0.75)
% yyaxis right
% plot(aqd.SSC2,'k'),hold on,plot(aqd.SSC1,'r')
% ylim([0 1000])

fid=fopen('Sept17_SmAqd.csv');
T=textscan(fid,'%*f %s %f %f','Delimiter',',','HeaderLines',2);fclose(fid);
hobo_cond=T{2};
hobo_temp=T{3};
hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');

aqd.hobo.cond=interp1(hobo_time,hobo_cond,aqd.time);
aqd.hobo.temp=interp1(hobo_time,hobo_temp,aqd.time);
aqd.hobo.sal=conduc2sali(aqd.hobo.cond/1000,aqd.hobo.temp,aqd.pres);

save('AQD_Sep17_SM','aqd')
%% Mar18 Main Channel
clear all,close all,clc
aqd=aqd_process('MMLG1801',0.05,2,2,[],[],0.5);
% offset for 11782 = 
offset=-9.35;

% remove bad data at beginning of deployment:
flds=fieldnames(aqd);
flds([11,19:21])=[];
for jj=1:length(flds)
    if jj<11 || jj>18
    aqd.(flds{jj})=aqd.(flds{jj})(1:end-1);
    else
    aqd.(flds{jj})=aqd.(flds{jj})(1:end-1,:);
    end
end

% remove air pressure signal using Met station and add HAB:
load('BogaleRiverWeather.mat')
Weather.AtmPres=fillmissing(Weather.AtmPres,'nearest');
airpres=interp1(Weather.datenum,Weather.AtmPres,aqd.time)./100;
aqd.depth=((aqd.pres-(airpres+offset)).*1.02)+0.10;
% aqd.depth(aqd.depth<3.5)=NaN;

% convert turbidity to SSC
aqd.SSC1=aqd.ext1.*0.122-34; % see obs calibration pdf
aqd.SSC2=aqd.ext2.*0.122-34;

fid=fopen('Mar18_LgChAqd.csv');
T=textscan(fid,'%*f %s %f %f','Delimiter',',','HeaderLines',2);fclose(fid);
hobo_cond=T{2};
hobo_temp=T{3};
hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');

aqd.hobo.cond=interp1(hobo_time,hobo_cond,aqd.time);
aqd.hobo.temp=interp1(hobo_time,hobo_temp,aqd.time);
aqd.hobo.sal=conduc2sali(aqd.hobo.cond/1000,aqd.hobo.temp,aqd.pres);

figure;
plot(aqd.depth)
figure;
plot(aqd.SSC1),hold on,plot(aqd.SSC2)
figure;
plot(aqd.hobo.sal)

%save('AQD_Mar18_LG','aqd')

%% Mar18 Small Channel
clear all,close all,clc
aqd=aqd_process('MMsm1801',0.05,2,2,[],[],0.5);
% offset for 6106 = 
offset=-9.91;

% remove bad data at beginning of deployment:
flds=fieldnames(aqd);
flds([11,19:21])=[];
for jj=1:length(flds)
    if jj<11 || jj>18
    aqd.(flds{jj})=aqd.(flds{jj})(46:1155);
    else
    aqd.(flds{jj})=aqd.(flds{jj})(46:1155,:);
    end
end

% remove air pressure signal using Met station and add HAB:
load('BogaleRiverWeather.mat')
Weather.AtmPres=fillmissing(Weather.AtmPres,'nearest');
airpres=interp1(Weather.datenum,Weather.AtmPres,aqd.time)./100;
aqd.depth=((aqd.pres-(airpres+offset)).*1.02)+0.10;
% aqd.depth(aqd.depth<3.5)=NaN;

% convert turbidity to SSC
aqd.SSC1=aqd.ext1.*0.111-24; % see obs calibration pdf
aqd.SSC2=aqd.ext2.*0.066-28;
aqd.SSC1(aqd.SSC1<20)=NaN;
aqd.SSC2(aqd.SSC2<20)=NaN;

fid=fopen('Mar18_SmChAqd.csv');
T=textscan(fid,'%*f %s %f %f','Delimiter',',','HeaderLines',2);fclose(fid);
hobo_cond=T{2};
hobo_temp=T{3};
hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');

aqd.hobo.cond=interp1(hobo_time,hobo_cond,aqd.time);
aqd.hobo.temp=interp1(hobo_time,hobo_temp,aqd.time);
aqd.hobo.sal=conduc2sali(aqd.hobo.cond/1000,aqd.hobo.temp,aqd.pres);
aqd.hobo.sal(aqd.hobo.sal<14)=NaN;

figure;
plot(aqd.depth)
figure;
plot(aqd.SSC1,'k'),hold on,plot(aqd.SSC2,'r')
figure;
plot(aqd.hobo.sal)

save('AQD_Mar18_SM','aqd')
%% Mar18 Agri
clear all,close all,clc
aqd=aqd_process('ag_ch01',0.05,2,2,[],[],0.5);
% offset for 6106 = 
offset=-9.35;

% remove bad data at beginning of deployment:
flds=fieldnames(aqd);
flds([11,19:21])=[];
for jj=1:length(flds)
    if jj<11 || jj>18
    aqd.(flds{jj})=aqd.(flds{jj})(7:582);
    else
    aqd.(flds{jj})=aqd.(flds{jj})(7:582,:);
    end
end

% remove air pressure signal using Met station and add HAB:
load('BogaleRiverWeather.mat')
Weather.AtmPres=fillmissing(Weather.AtmPres,'nearest');
airpres=interp1(Weather.datenum,Weather.AtmPres,aqd.time)./100;
aqd.depth=((aqd.pres-(airpres+offset)).*1.02)+0.10;

% convert turbidity to SSC
aqd.SSC1=aqd.ext1.*0.122-34; % see obs calibration pdf
aqd.SSC2=aqd.ext2.*0.122-34;
aqd.SSC1(aqd.depth<0.24 | aqd.SSC1<20)=NaN;
aqd.SSC2(aqd.depth<0.54| aqd.SSC2<20)=NaN;

fid=fopen('Mar18_AgFieldAqd.csv');
T=textscan(fid,'%*f %s %f %f','Delimiter',',','HeaderLines',2);fclose(fid);
hobo_cond=T{2};
hobo_temp=T{3};
hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');

aqd.hobo.cond=interp1(hobo_time,hobo_cond,aqd.time);
aqd.hobo.temp=interp1(hobo_time,hobo_temp,aqd.time);
aqd.hobo.sal=conduc2sali(aqd.hobo.cond/1000,aqd.hobo.temp,aqd.pres);
aqd.hobo.sal(aqd.depth<0.24)=NaN;

aqd.depth(aqd.depth<0.1)=NaN;


figure;
plot(aqd.depth)
yyaxis right
plot(aqd.SSC1,'k'),hold on,plot(aqd.SSC2,'r')
figure;
plot(aqd.depth)
yyaxis right
plot(aqd.hobo.sal)

save('AQD_Mar18_AG','aqd')
%% Sept19 Main Ch
clear all,close all,clc
aqd=aqd_process('BigCh01',0.05,2,2,[],[],0.5);
% offset for 11782 = 
offset=-9.35;

% remove bad data at beginning of deployment:
flds=fieldnames(aqd);
% flds([11,19:21])=[];
% for jj=1:length(flds)
%     if jj<11 || jj>18
%     aqd.(flds{jj})=aqd.(flds{jj})(96:3073);
%     else
%     aqd.(flds{jj})=aqd.(flds{jj})(96:3073,:);
%     end
% end

% remove air pressure signal using Met station and add HAB:
load('BogaleRiverWeather.mat')
Weather.AtmPres=fillmissing(Weather.AtmPres,'nearest');
airpres=interp1(Weather.datenum,Weather.AtmPres,aqd.time)./100;
aqd.depth=((aqd.pres-(airpres+offset)).*1.02)+0.10;

% convert turbidity to SSC
aqd.SSC1=aqd.ext1.*0.067-34; % see obs calibration pdf
aqd.SSC2=aqd.ext2.*0.067-34;
ssc_diff=abs([aqd.SSC2(2:end);0]-aqd.SSC2);
aqd.SSC2(ssc_diff>200)=NaN;
aqd.SSC1(aqd.SSC1<20 | aqd.SSC1>1500)=NaN;
aqd.SSC2(aqd.SSC2<20| aqd.SSC2>1500)=NaN;

%FIND CONDUCTIVITY DATA
% fid=fopen('Mar18_AgFieldAqd.csv');
% T=textscan(fid,'%*f %s %f %f','Delimiter',',','HeaderLines',2);fclose(fid);
% hobo_cond=T{2};
% hobo_temp=T{3};
% hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');
% 
% aqd.hobo.cond=interp1(hobo_time,hobo_cond,aqd.time);
% aqd.hobo.temp=interp1(hobo_time,hobo_temp,aqd.time);
% aqd.hobo.sal=conduc2sali(aqd.hobo.cond/1000,aqd.hobo.temp,aqd.pres);
% aqd.hobo.sal(aqd.depth<0.24)=NaN;

figure;
plot(aqd.depth)
yyaxis right
plot(aqd.SSC1,'k'),hold on,plot(aqd.SSC2,'r')

save('AQD_Sept19_LG','aqd')

%% Sept19 Ag Ch
clear all,close all,clc
aqd=aqd_process('AgCh1901',0.05,2,2,[],[],0.5);
% offset for 6106 = 
offset=-9.91;

% remove bad data at beginning of deployment:
flds=fieldnames(aqd);
flds([11,19:21])=[];
for jj=1:length(flds)
    if jj<11 || jj>18
    aqd.(flds{jj})=aqd.(flds{jj})(289:4096);
    else
    aqd.(flds{jj})=aqd.(flds{jj})(289:4096,:);
    end
end

% remove air pressure signal using Met station and add HAB:
load('BogaleRiverWeather.mat')
Weather.AtmPres=fillmissing(Weather.AtmPres,'nearest');
airpres=interp1(Weather.datenum,Weather.AtmPres,aqd.time)./100;
aqd.depth=((aqd.pres-(airpres+offset)).*1.02)+0.10;

% convert turbidity to SSC
aqd.SSC1=aqd.ext1.*0.111-24; % see obs calibration pdf
aqd.SSC2=aqd.ext2.*0.066-28;
aqd.SSC1(aqd.depth<0.33 | aqd.SSC1>1000)=NaN;
aqd.SSC2(aqd.depth<0.53)=NaN;

% load conductivity
fid=fopen('AgriCh_092719.csv');
T=textscan(fid,'%*f %s %f %f','Delimiter',',','HeaderLines',2);fclose(fid);
hobo_cond=T{2};
hobo_temp=T{3};
hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');

aqd.hobo.cond=interp1(hobo_time,hobo_cond,aqd.time);
aqd.hobo.temp=interp1(hobo_time,hobo_temp,aqd.time);
aqd.hobo.temp(aqd.hobo.temp>35)=35;
aqd.hobo.sal=conduc2sali(aqd.hobo.cond/1000,aqd.hobo.temp,aqd.pres);
aqd.hobo.sal(aqd.depth<0.33)=NaN;
clear T hobo*

% load pressure in field:
fid=fopen('AgriField_Pressure_092719.csv');
T=textscan(fid,'%*f %s %f %f','Delimiter',',','HeaderLines',2);fclose(fid);
hobo_pres=T{2}; 
hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');
aqd.field.pres=interp1(hobo_time,hobo_pres,aqd.time)./100;% convert from mbar to dbar
aqd.field.depth=aqd.field.pres-airpres;
aqd.field.depth(aqd.field.depth<0.03)=NaN;

aqd.depth(aqd.depth<0.1)=NaN;

figure;
plot(aqd.depth)
yyaxis right
plot(aqd.SSC1,'k'),hold on,plot(aqd.SSC2,'r')
figure;
plot(aqd.depth)
yyaxis right
plot(aqd.hobo.sal)

save('AQD_Sept19_AG','aqd')
%% Sept19 Sm CH (high freq aqd...SN:11987)
clear all,close all,clc
aqd=aqd_process('SMCH1916',0.05,1,2,[],[],0.5);
% offset for 11987 = 
offset=-9.6;

% remove bad data at beginning of deployment:
flds=fieldnames(aqd);
flds([20:21])=[];
for jj=1:length(flds)
    [~,nn]=size(aqd.(flds{jj}));
    if nn==1
    aqd.(flds{jj})=aqd.(flds{jj})(53:598);
    else
    aqd.(flds{jj})=aqd.(flds{jj})(53:598,:);
    end
end

% remove air pressure signal using Met station and add HAB:
load('BogaleRiverWeather.mat')
Weather.AtmPres=fillmissing(Weather.AtmPres,'nearest');
airpres=interp1(Weather.datenum,Weather.AtmPres,aqd.time)./100;
aqd.depth=((aqd.pres-(airpres+offset)).*1.02)+0.10;

% convert turbidity to SSC (use approx for this cal)
aqd.SSC1=NaN(length(aqd.ext1),1);
aqd.SSC2=aqd.ext2.*0.07-43; 
aqd.SSC2(aqd.depth<0.33)=NaN;

figure;
plot(aqd.depth),
yyaxis right
plot(aqd.SSC2)

save('AQD_Sep19_SM','aqd')