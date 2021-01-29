clear all,close all,clc

% load aqd data output from the Nortek conversion software
cd C:\GLOVER\data\myanmar\AQD\BR_Aqd_Sept19
aqd=aqdlr2x('AgCh1901.hdr');
aqd.header.transmat = [1.5774 -0.7891 -0.7891;...
    0.0000 -1.3662 1.3662;...
    0.3677 0.3677 0.3677];

% add met station data:
load('C:\GLOVER\output\myanmar\longterminst\BogaleRiverWeather.mat')
[~,ia,iw]=intersect(aqd.time,Weather.datenum);
aqd.atm_pres = NaN(size(aqd.time));
aqd.atm_pres(ia) = Weather.AtmPres(iw)./100;
aqd.atm_pres=fillmissing(aqd.atm_pres,'nearest');

% figure;plot(aqd.pres),hold on,plot(aqd.atm_pres)
% find the pressure offset between aqd and met station
offset = nanmean(aqd.atm_pres(1:90)) - nanmean(aqd.pres(1:90));
aqd.depth = aqd.pres+offset-aqd.atm_pres;
aqd.depth(aqd.depth<0.25)=NaN;

% add instrument height to get total water depth:
aqd.depth = aqd.depth+0.15;
%add sensor head elevation to make z=MAB
aqd.z=aqd.z+0.15; 

%
fields=fieldnames(aqd);
%remove bad data at beginning and end, and nan some bad data intervals
for jj=[1:10,21:22]
    % remove data from 1D data fields
    aqd.(fields{jj})([1:292,4096:end])=[];

% figure; plot(aqd.(fields{jj})),title(fields{jj})
end
for jj=13:18
    % remove data from 2D data fields
    aqd.(fields{jj})([1:292,4096:end],:)=[];
end
% calculate water slope in m/hr
deltaT=round((1/24)./(aqd.time(2)-aqd.time(1)));
aqd.slope=(aqd.depth(2:end)-aqd.depth(1:end-1)).*deltaT;
aqd.slope(end+1)=aqd.slope(end);
aqd.slope(abs(aqd.slope)>1.2)=NaN;
aqd.slope = movmean(aqd.slope,11);

% calculate SSC from the OBS counts
% convert counts to mV
aqd.ext1=(aqd.ext1/65536)*5*1000;
aqd.ext2=(aqd.ext2/65536)*5*1000;
% then convert mV to SSC;AGD6106
aqd.ssc1 = aqd.ext1*0.11;
aqd.ssc2 = aqd.ext2*0.066;
aqd.ssc1(aqd.ssc1<2)=NaN;aqd.ssc2(aqd.ssc2<2)=NaN;
aqd.ssc1(3500:end)=NaN;
aqd.ssc1=movmean(aqd.ssc1,3);aqd.ssc2=movmean(aqd.ssc2,3);

figure;plot(aqd.ssc1,'k'),hold on,plot(aqd.ssc2,'r')
%%
% rotate the velocity components to ENU on head/roll/tilt
for jj=1:length(aqd.time)
    [aqd.v1(jj,:),aqd.v2(jj,:),aqd.v3(jj,:)]=xyz2enu(...
        aqd.v1(jj,:),aqd.v2(jj,:),aqd.v3(jj,:),...
        aqd.head(jj),aqd.pitch(jj),aqd.roll(jj),aqd.header.transmat,0,0);
end

% find out of water values:
aqd.air=repmat(aqd.z,length(aqd.time),1)>repmat(aqd.depth,1,aqd.header.numCells);
aqd.air(isnan(aqd.depth),:)=1;
for jj=13:18
    aqd.(fields{jj})(aqd.air)=NaN;
end

% depth avg the velocity
aqd.v1_mean=nanmean(aqd.v1,2);
aqd.v2_mean=nanmean(aqd.v2,2);
aqd.v3_mean=nanmean(aqd.v3,2);

% calculate enu speed and direction
[aqd.spd,aqd.dir]=uv2sd(aqd.v1,aqd.v2,aqd.v3);

% calculate mean speed and dir:
aqd.spd_mean=nanmean(aqd.spd,2);
aqd.dir_mean=nanmean(aqd.dir,2);

% add salinity data:
%add conductivity
fid=fopen('AgriCh_092719.csv');
T=textscan(fid,'%*f %s %f %*f','Delimiter',',','HeaderLines',2);
fclose(fid);
% convert conductivity to salinity
hobo_sal=T{2};
hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');

[~,ih,ia]=intersect(hobo_time,aqd.time);
aqd.sal=NaN([length(aqd.time),1]);
aqd.sal(ia)=hobo_sal(ih);aqd.sal=fillmissing(aqd.sal,'nearest');
aqd.sal=conduc2sali(aqd.sal/1000,30*ones(size(aqd.temp)),aqd.depth);

%%
figure;
subplot(311),pcolor(aqd.time,aqd.z,aqd.v1'),shading flat,colorbar
yyaxis right,plot(aqd.time,aqd.roll)
subplot(312),pcolor(aqd.time,aqd.z,aqd.v2'),shading flat,colorbar
yyaxis right,plot(aqd.time,aqd.roll)
subplot(313),pcolor(aqd.time,aqd.z,aqd.v3'),shading flat,colorbar
yyaxis right,plot(aqd.time,aqd.roll)

figure;
subplot(211),pcolor(aqd.time,aqd.z,aqd.spd'),shading flat,colorbar,caxis([0 1])
subplot(212),pcolor(aqd.time,aqd.z,aqd.dir'),shading flat,colorbar

figure;
plot(aqd.dir(:,2),'k'),hold on
plot(aqd.dir(:,20),'r')
yyaxis right
plot(aqd.depth)
figure;
subplot(311),plot(aqd.head),ylabel('head')
subplot(312),plot(aqd.pitch),ylabel('pitch')
subplot(313),plot(aqd.roll),ylabel('roll')


%% go to output data folder

cd C:\GLOVER\output\myanmar\aqd
save('aqd_sep19_ag','aqd')
