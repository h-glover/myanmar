clear all,close all,clc

% load aqd data output from the Nortek conversion software
cd C:\GLOVER\data\myanmar\AQD\BR_Aqd_Sept19
load('SMCH1916_burstdata.mat')

% Make the time vector in the aqd and remove the extraneous time fields
F=fields(fullaqd);
for jj=1:19
    aqd.(F{jj}) = double(horzcat(fullaqd.(F{jj})));
    
    if jj<7
        aqd.(F{jj}) = aqd.(F{jj})(1,:)';
    end
end
aqd.time=datenum(aqd.year,aqd.month,aqd.day,aqd.hour,aqd.minute,aqd.second);
aqd=rmfield(aqd,{'month','day','year','hour','minute','second'});

% fix the turbidity
aqd.ext1([1,512],:)=NaN;
aqd.ext2([1,512],:)=NaN;
aqd.ext1=(aqd.ext1/65536)*5*1000;
aqd.ext2=(aqd.ext2/65536)*5*1000;
% then convert mV to SSC;
aqd.ssc1 = aqd.ext1*0.06;
aqd.ssc1 = nanmean(aqd.ssc1)';
%d = abs((aqd.ssc1(3:end)-aqd.ssc1(1:end-2))/2);d=[0;d;0];
aqd.ssc1(aqd.ssc1<25 | aqd.ssc1>270)=NaN;
aqd.ssc2 = aqd.ext2*0.06;
% d = abs((aqd.ssc2(3:end)-aqd.ssc2(1:end-2))/2);d=[0;d;0];
aqd.ssc2(aqd.ssc2<25 | aqd.ssc2>270)=NaN;
aqd.ssc2 = nanmean(aqd.ssc2)';

% average all values that aren't velocity
for jj=7:19
    aqd.(F{jj}) = nanmean(aqd.(F{jj}));
    aqd.(F{jj})=aqd.(F{jj})';
end


% add met station data:
load('C:\GLOVER\output\myanmar\longterminst\BogaleRiverWeather.mat')
[~,ia,iw]=intersect(aqd.time,Weather.datenum);
aqd.atm_pres = NaN(size(aqd.time));
aqd.atm_pres(ia) = Weather.AtmPres(iw)./100;
aqd.atm_pres=fillmissing(aqd.atm_pres,'nearest');

% figure;plot(aqd.pres),hold on,plot(aqd.atm_pres)
% find the pressure offset between aqd and met station
offset = nanmean(aqd.atm_pres(1:50)) - nanmean(aqd.pres(1:50));
aqd.depth = aqd.pres+offset-aqd.atm_pres;
aqd.depth(aqd.depth<0.25)=NaN;

% add instrument height to get total water depth:
aqd.depth = aqd.depth+0.15;
%add sensor head elevation to make z=MAB
% aqd.z=aqd.z+0.15; 

% remove bad data at beginning and end, and nan some bad data intervals
fields=fieldnames(aqd);
for jj=1:length(fields)
    % remove data from 1D data fields
    aqd.(fields{jj})([1:53,598:end])=[];

% figure; plot(aqd.(fields{jj})),title(fields{jj})
end
fullaqd=rmfield(fullaqd,F([1:12,24:28]));
fullaqd([1:53,598:end])=[];

% calculate water slope in m/hr
deltaT=round((1/24)./(aqd.time(2)-aqd.time(1)));
aqd.slope=(aqd.depth(2:end)-aqd.depth(1:end-1)).*deltaT;
aqd.slope(end+1)=aqd.slope(end);
aqd.slope(abs(aqd.slope)>1.2)=NaN;
aqd.slope = movmean(aqd.slope,11);


%
% make a new vector for the burst-organized velocity data
% row = within burst, col = depth cell, 3=time
aqd.v1b=double(cat(3,fullaqd.v1));
aqd.v2b=double(cat(3,fullaqd.v2));
aqd.v3b=double(cat(3,fullaqd.v3));

% reshape to row = time step, col = within burst, 3=depth bin
aqd.v1b = permute(aqd.v1b,[3,1,2]);
aqd.v2b = permute(aqd.v2b,[3,1,2]);
aqd.v3b = permute(aqd.v3b,[3,1,2]);

% remove bad vals from instrument movement
aqd.v1b([63:115, 205:240,499:514],:,:)=NaN;
aqd.v2b([63:115, 205:240,499:514],:,:)=NaN;
aqd.v3b([63:115, 205:240,499:514],:,:)=NaN;

% rotate the velocity components to along,across,up based on tilt and
% channel orientation
% probes were all oriented upwards; pos z is upward, pos u is northward
% coord sys = XYZ where instrument was pointed N
% first conveert xyz to beam to enu
aqd.header.transmat = [1.5774 -0.7891 -0.7891;...
    0.0000 -1.3662 1.3662;...
    0.3677 0.3677 0.3677];
% enu conversion for within burst velocity
% use burst avg head, pitch, roll and compute for each burst
[L,aqd.header.spb,aqd.header.bins] =size(aqd.v1b); 
for kk=1:aqd.header.bins
for jj=1:L
    [aqd.v1b(jj,:,kk),aqd.v2b(jj,:,kk),aqd.v3b(jj,:,kk)]=xyz2enu(...
        aqd.v1b(jj,:,kk),aqd.v2b(jj,:,kk),aqd.v3b(jj,:,kk),...
        aqd.heading(jj),aqd.pitch(jj),aqd.roll(jj),aqd.header.transmat/4096,0,0);
end
end

% % despike the burst velocity data
% % [VelX{i}(:,j),ip] = func_despike_phasespace3d(VelX{i}(:,j),0,2);
% for kk=1:aqd.header.bins
% for jj=1:L
%     [aqd.v1b(jj,:,kk), ~] = func_despike_phasespace3d(aqd.v1b(jj,:,kk), 0, 0);
%     [aqd.v2b(jj,:,kk), ~] = func_despike_phasespace3d(aqd.v2b(jj,:,kk), 0, 0);
%     [aqd.v3b(jj,:,kk), ~] = func_despike_phasespace3d(aqd.v3b(jj,:,kk), 0, 0);
% end
% end

% burst avg the velocity
aqd.v1=permute(nanmean(aqd.v1b,2),[1,3,2]);
aqd.v2=permute(nanmean(aqd.v2b,2),[1,3,2]);
aqd.v3=permute(nanmean(aqd.v3b,2),[1,3,2]);

% find out of water values:
% nowacki: edited so aqd.z represents midpoint of cell, not bottom.
aqd.header.blank = 0.055; %m blanking distance
aqd.header.dz = 0.05; %m bin height
aqd.z = 0.15 + (aqd.header.blank + (0:(aqd.header.bins-1)) * aqd.header.dz + aqd.header.dz/2); 

aqd.air=repmat(aqd.z,length(aqd.time),1)>repmat(aqd.depth,1,aqd.header.bins);
aqd.air(isnan(aqd.depth),:)=1;
aqd.v1(aqd.air)=NaN;
aqd.v2(aqd.air)=NaN;
aqd.v3(aqd.air)=NaN;

% calculate enu speed and direction
[aqd.spd,aqd.dir]=uv2sd(aqd.v1,aqd.v2,aqd.v3);

% depth avg the velocity
aqd.v1_mean=nanmean(aqd.v1,2);
aqd.v2_mean=nanmean(aqd.v2,2);
aqd.v3_mean=nanmean(aqd.v3,2);

% calculate mean speed and dir:
aqd.spd_mean=nanmean(aqd.spd,2);
aqd.dir_mean=nanmean(aqd.dir,2);

% % add salinity data:
% %add conductivity
% fid=fopen('AgriCh_092719.csv');
% T=textscan(fid,'%*f %s %f %*f','Delimiter',',','HeaderLines',2);
% fclose(fid);
% % convert conductivity to salinity
% hobo_sal=T{2};
% hobo_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');
% 
% [~,ih,ia]=intersect(hobo_time,aqd.time);
% aqd.sal=NaN([length(aqd.time),1]);
% aqd.sal(ia)=hobo_sal(ih);aqd.sal=fillmissing(aqd.sal,'nearest');
% aqd.sal=conduc2sali(aqd.sal/1000,30*ones(size(aqd.temp)),aqd.depth);


%%
figure;
subplot(311),pcolor(aqd.time,aqd.z,aqd.v1'),shading flat,colorbar
hold on,plot(aqd.time,aqd.depth,'k')
subplot(312),pcolor(aqd.time,aqd.z,aqd.v2'),shading flat,colorbar
hold on,plot(aqd.time,aqd.depth,'k')
subplot(313),pcolor(aqd.time,aqd.z,aqd.v3'),shading flat,colorbar
hold on,plot(aqd.time,aqd.depth,'k')

figure;
subplot(211),pcolor(aqd.time,aqd.z,aqd.spd'),shading flat,colorbar,caxis([0 1])
subplot(212),pcolor(aqd.time,aqd.z,aqd.dir'),shading flat,colorbar

figure;
subplot(311),plot(aqd.time,aqd.spd_mean,'r')
yyaxis right,plot(aqd.time,aqd.depth,'k')
subplot(312),plot(aqd.time,aqd.dir_mean,'r')
yyaxis right,plot(aqd.time,aqd.depth,'k')
subplot(313),quiver(aqd.time,zeros(size(aqd.time)),aqd.v1_mean,aqd.v2_mean)

figure;
subplot(311),plot(aqd.heading),ylabel('head')
subplot(312),plot(aqd.pitch),ylabel('pitch')
subplot(313),plot(aqd.roll),ylabel('roll')


%% go to output data folder

cd C:\GLOVER\output\myanmar\aqd
save('aqd_sep19_lc','aqd')
