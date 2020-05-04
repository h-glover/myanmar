% Process aquadopp data from MeinMaHla Island
clear all,close all,clc

aqd=aqd_process('ag_ch01',0.05,2,2,[],[],0.5);
aqd.time=datevec(aqd.time);aqd.time=datenum(aqd.time);
aqd.pres=aqd.pres-0.5; %0.47, 0.7
% STOP and use the code from aqdpressuretrim to cut out first 7 and last 3 values
% ext1=20cmab, ext2=50cmab, salinity=20cmab
aqd.ext1(aqd.pres<0.20)=NaN;
aqd.ext2(aqd.pres<0.50)=NaN;

fid=fopen('Mar18_AgFieldAqd.csv');
T=textscan(fid,'%*f %s %f %f','Delimiter',',','HeaderLines',2);fclose(fid);
hobo_sal=T{2};
hobo_temp=T{3};

hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');
hobo_time=datevec(hobo_time);
hobo_time(:,5)=round(hobo_time(:,5)/5)*5;
hobo_time=datenum(hobo_time);

hobo_5min=unique(hobo_time);
for jj=1:length(hobo_5min)
    hobo_5sal(jj)=nanmean(hobo_sal(hobo_time==hobo_5min(jj)));
    hobo_5temp(jj)=nanmean(hobo_temp(hobo_time==hobo_5min(jj)));
end

[~,~,ih]=intersect(aqd.time,hobo_5min);
aqd.hobo_sal=hobo_5sal(ih)';
aqd.hobo_temp=hobo_5temp(ih)';
aqd.hobo_sal([1,end])=NaN;
aqd.hobo_temp([1,end])=NaN;
aqd.sal=conduc2sali(aqd.hobo_sal/1000,aqd.temp,aqd.pres);
aqd.sal(aqd.pres<0.50)=NaN;

save('MMAG_Mar18_aqd','aqd')

%%

clear all,close all,clc

aqd=aqd_process('MMsm1801',0.05,2,2,[],[],0.5);
aqd.time=datevec(aqd.time);aqd.time=datenum(aqd.time);
aqd.pres=aqd.pres-0.2;
plot(aqd.pres)
aqd=aqdpressuretrim(aqd,0.2);
% ext1=20cmab, ext2=50cmab, salinity=20cmab
aqd.ext1(aqd.pres<0.20)=NaN;
aqd.ext2(aqd.pres<0.50)=NaN;

fid=fopen('Mar18_SmChAqd.csv');
T=textscan(fid,'%*f %s %f %f','Delimiter',',','HeaderLines',2);fclose(fid);
hobo_sal=T{2};
hobo_temp=T{3};

hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');
hobo_time=datevec(hobo_time);
hobo_time(:,5)=round(hobo_time(:,5)/5)*5;
hobo_time=datenum(hobo_time);

hobo_5min=unique(hobo_time);
for jj=1:length(hobo_5min)
    hobo_5sal(jj)=nanmean(hobo_sal(hobo_time==hobo_5min(jj)));
    hobo_5temp(jj)=nanmean(hobo_temp(hobo_time==hobo_5min(jj)));
end

[~,~,ih]=intersect(aqd.time,hobo_5min);
aqd.hobo_sal=hobo_5sal(ih)';
aqd.hobo_temp=hobo_5temp(ih)';
aqd.hobo_sal([1,end])=NaN;
aqd.hobo_temp([1,end])=NaN;
aqd.sal=conduc2sali(aqd.hobo_sal/1000,aqd.temp,aqd.pres);
aqd.sal(aqd.pres<0.50)=NaN;

save('MMsm_Mar18_aqd','aqd')
%%
clear all,close all,clc

aqd=aqd_process('MMLG1801',0.05,2,2,[],[],0.5);
aqd.time=datevec(aqd.time);aqd.time=datenum(aqd.time);
plot(aqd.pres)
aqd.pres=aqd.pres-0.5; %same offset as ag ch deployment
aqd=aqdpressuretrim(aqd,2);

fid=fopen('Mar18_LgChAqd.csv');
T=textscan(fid,'%*f %s %f %f','Delimiter',',','HeaderLines',2);fclose(fid);
hobo_sal=T{2};
hobo_temp=T{3};

hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');
hobo_time=datevec(hobo_time);
hobo_time(:,5)=round(hobo_time(:,5)/5)*5;
hobo_time=datenum(hobo_time);

hobo_5min=unique(hobo_time);
for jj=1:length(hobo_5min)
    hobo_5sal(jj)=nanmean(hobo_sal(hobo_time==hobo_5min(jj)));
    hobo_5temp(jj)=nanmean(hobo_temp(hobo_time==hobo_5min(jj)));
end

[~,~,ih]=intersect(aqd.time,hobo_5min);
aqd.hobo_sal=hobo_5sal(ih)';
aqd.hobo_temp=hobo_5temp(ih)';
aqd.hobo_sal([1,end])=NaN;
aqd.hobo_temp([1,end])=NaN;
aqd.sal=conduc2sali(aqd.hobo_sal/1000,aqd.temp,aqd.pres);

save('MMLG_Mar18_aqd','aqd')
% %%
% clear all,close all,clc
% load('MMLG_Mar18_aqd')
% 
% 
% 
% figure;
% subplot(221)
% scatter(aqd.ext1,aqd.sal),title('ext1 v sal')
% subplot(222)
% scatter(aqd.ext1,aqd.temp),title('ext1 v temp')
% subplot(223)
% scatter(aqd.ext2,aqd.sal),title('ext2 v sal')
% subplot(224)
% scatter(aqd.ext2,aqd.temp),title('ext2 v temp')