%% load enx data into adcp-type structure:

clear all,close all,clc
cd C:\GLOVER\data\myanmar\transecdischarge\YR_Spring_Jan19\raw_adcp


F=dir('ENX*.mat');
for ff=1:length(F)
    load(F(ff).name)
    
    adcp(ff).binsize=RDIBinSize;
    adcp(ff).bindepth=RDIBin1Mid:RDIBinSize:(RDIBinSize*100+RDIBin1Mid-RDIBinSize);

    % make timevec
    adcp(ff).time=datenum(2000+SerYear,SerMon,SerDay,SerHour,SerMin,SerSec);
    
    % add in boat lat and lon
    adcp(ff).lat=AnFLatDeg;
    adcp(ff).lon=AnFLonDeg;
    
    % add in boat E and N vectors:
    adcp(ff).boat_east=AnNVEmmpersec/10;
    adcp(ff).boat_north=AnNVNmmpersec/10;
    
    % convert vel to cm/s
    adcp(ff).east=(SerEmmpersec/10); % beam 1
    adcp(ff).north=(SerNmmpersec/10); % beam 2
    adcp(ff).up=SerVmmpersec/10; % beam 3
    % remove low corr values and erroneous values
    SerC1cnt=100*SerC1cnt./255;
    SerC2cnt=100*SerC2cnt./255;
    SerC3cnt=100*SerC3cnt./255;
    adcp(ff).east(SerC1cnt<30)=NaN;
    adcp(ff).north(SerC2cnt<30)=NaN;
    adcp(ff).up(SerC3cnt<30)=NaN;
    
    adcp(ff).east(abs(adcp(ff).east)>500)=NaN;
    adcp(ff).north(abs(adcp(ff).north)>500)=NaN;
    adcp(ff).up(abs(adcp(ff).up)>500)=NaN;
    
    
    % remove everything below 22 m depth
    adcp(ff).east(:,adcp(ff).bindepth>22)=[];
    adcp(ff).north(:,adcp(ff).bindepth>22)=[];
    adcp(ff).up(:,adcp(ff).bindepth>22)=[];
    adcp(ff).bindepth(adcp(ff).bindepth>22)=[];
    
    % rotate to match other adcp data:
    adcp(ff).east=adcp(ff).east';
    adcp(ff).north=adcp(ff).north';
    adcp(ff).up=adcp(ff).up';
    
    t_srt(ff)=adcp(ff).time(1);
 
end
[~,srt]=sort(t_srt,'ascend');
adcp=adcp(srt);

cd C:\GLOVER\output\myanmar\transectdischarge\yr_spring_jan19
save('YR_spring_jan19','adcp')
