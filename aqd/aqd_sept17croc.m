clear all,close all,clc

% load aqd data output from the Nortek conversion software
cd C:\GLOVER\data\myanmar\AQD\BR_Aqd_Sept17\MainCrocStation
aqd=aqdlr2x('BRmain01.hdr');
aqd.header.transmat = [1.5774 -0.7891 -0.7891;...
    0.0000 -1.3662 1.3662;...
    0.3677 0.3677 0.3677];

%add rbr data
load('065625_20170913_1535.mat')
RBR.time = datenum(RBR.sampletimes,'yyyy-mm-dd HH:MM:SS.FFF');

% add met station data:
load('C:\GLOVER\output\myanmar\longterminst\BogaleRiverWeather.mat')
% find the pressure offset between aqd and met station and RBR

[~,ia,ir]=intersect(aqd.time,RBR.time);
aqd.rbr_pres = NaN(size(aqd.time));
aqd.rbr_pres(ia) = RBR.data(ir,3);
aqd.rbr_pres=fillmissing(aqd.rbr_pres,'linear');
[~,ia,iw]=intersect(aqd.time,Weather.datenum);
aqd.atm_pres = NaN(size(aqd.time));
aqd.atm_pres(ia) = Weather.AtmPres(iw)./100;
aqd.atm_pres=fillmissing(aqd.atm_pres,'linear');

% fix the pressure using the atmospheric pressure and the RBR record (offset in sensor)
idx = 4633:4683; %idx1 = 1:70;
aqd.depth=aqd.pres + nanmean(aqd.rbr_pres(idx)-aqd.pres(idx))- aqd.atm_pres;
% add instrument height to get total water depth:
aqd.depth = aqd.depth+0.15;
%add sensor head elevation to make z=MAB
aqd.z=aqd.z+0.15; 

%%
fields=fieldnames(aqd);
%remove bad data at beginning and end, and nan some bad data intervals
for jj=[1:10,21:23]
    % remove data from 1D data fields
    aqd.(fields{jj})([1:480,4625:4683])=[];
    if jj>1
        aqd.(fields{jj})([569:765,3665:3787])=NaN;
    end
end
for jj=13:18
    % remove data from 2D data fields
    aqd.(fields{jj})([1:480,4625:4683],:)=[];
    aqd.(fields{jj})([569:765,3665:3787],:)=NaN;
end

% nan some bad data intervals
% % hand selected bad values due to mud? on sensor
% aqd.depth([2180:2300,2925:3050,3690:3800,4412:4550,5155:5290,...
%     5892:6036,6627:6761])=NaN;

% calculate water slope in m/hr
deltaT=round((1/24)./(aqd.time(2)-aqd.time(1)));
aqd.slope=(aqd.depth(2:end)-aqd.depth(1:end-1)).*deltaT;
aqd.slope(end+1)=aqd.slope(end);
aqd.slope(abs(aqd.slope)>1.2)=NaN;
aqd.slope = movmean(aqd.slope,11);
%% calculate SSC from the OBS counts
% convert counts to mV
% aqd.ext1(aqd.depth<0.17 | isnan(aqd.depth))=NaN;
aqd.ext1=(aqd.ext1/65536)*5*1000;
aqd.ext2=(aqd.ext2/65536)*5*1000;
% then convert mV to SSC; ext1 = 9192 = cal(2)
aqd.ssc1 = aqd.ext1*0.12-1;
aqd.ssc2 = aqd.ext2*0.12-1;

aqd.ssc1(aqd.ssc1<1)=NaN;
aqd.ssc2(aqd.ssc2<1)=NaN;
aqd.ssc1(abs(aqd.ssc1(1:end-1)-aqd.ssc1(2:end))>70)=NaN;
aqd.ssc2(abs(aqd.ssc2(1:end-1)-aqd.ssc2(2:end))>70)=NaN;


%% rotate the velocity components to ENU on head/roll/tilt
% for jj=1:length(aqd.time)
%     [aqd.v1(jj,:),aqd.v2(jj,:),aqd.v3(jj,:)]=xyz2enu(...
%         aqd.v1(jj,:),aqd.v2(jj,:),aqd.v3(jj,:),...
%         aqd.head(jj),aqd.pitch(jj),aqd.roll(jj),aqd.header.transmat,0,0);
% end

% find out of water values:
aqd.depth(aqd.depth<=0.5)=NaN;
aqd.air=repmat(aqd.z,length(aqd.time),1)>repmat(aqd.depth,1,aqd.header.numCells);
aqd.air(isnan(aqd.depth),:)=1;
for jj=13:18
    aqd.(fields{jj})(aqd.air)=NaN;
end

% depth avg the velocity
aqd.v1_mean=nanmean(aqd.v1,2);
aqd.v2_mean=nanmean(aqd.v2,2);
aqd.v3_mean=nanmean(aqd.v3,2);
% figure;
% quiver(aqd.time,zeros([1,length(aqd.time)]),aqd.v1_mean,aqd.v2_mean)
% yyaxis right
% plot(aqd.time,aqd.depth)

% calculate enu speed and direction
[aqd.spd,aqd.dir]=uv2sd(aqd.v1,aqd.v2,aqd.v3);

% calculate mean speed and dir:
aqd.spd_mean=nanmean(aqd.spd,2);
aqd.dir_mean=nanmean(aqd.dir,2);

%% add other RBR data??

%%
% figure;
% subplot(311),pcolor(aqd.time,aqd.z,aqd.v1'),shading flat,colorbar
% caxis([-1.5 1.5]),title('Croc, Sept 2017'),ylabel('v1')
% subplot(312),pcolor(aqd.time,aqd.z,aqd.v2'),shading flat,colorbar
% caxis([-1.5 1.5]),ylabel('v2')
% subplot(313),pcolor(aqd.time,aqd.z,aqd.v3'),shading flat,colorbar
% caxis([-1.5 1.5]),ylabel('v3')
% 
% figure;
% subplot(211),pcolor(aqd.time,aqd.z,aqd.spd'),shading flat,colorbar
% caxis([0 2.5]),ylabel('spd'),title('Croc, Sept 2017')
% subplot(212),pcolor(aqd.time,aqd.z,aqd.dir'),shading flat,colorbar
% caxis([0 360]),ylabel('dir')
% % 
% % figure;
% % plot(aqd.dir(:,2),'k'),hold on
% % plot(aqd.dir(:,20),'r')
% % yyaxis right
% % plot(aqd.depth)
% figure;
% subplot(311),plot(aqd.head),ylabel('head')
% subplot(312),plot(aqd.pitch),ylabel('pitch')
% subplot(313),plot(aqd.roll),ylabel('roll')

%% go to output data folder

cd C:\GLOVER\output\myanmar\aqd
save('aqd_sep17_hc','aqd')

