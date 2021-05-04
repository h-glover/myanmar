%% make a struc that is just from interior island surveys
clear all,close all,clc

cd C:\GLOVER\data\myanmar\mmi_velocity\BR_Int

% Folder containing only text files of each transect
FileList=dir('*ASC.TXT');

for n = 1:length(FileList)
    disp(['processing ' num2str(n)]);
    
    adcp(n) = adcp_ascii(FileList(n).name);

end

for jj=1:length(adcp)
        adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
        adcp(jj).hour,adcp(jj).minute,adcp(jj).second);
end
cd C:\GLOVER\output\myanmar\mmi_discharge
save('br_08mar2018','adcp')
%%
clear all,close all,clc

cd C:\GLOVER\data\myanmar\mmi_velocity\EastChannel_Mar18
% Folder containing only text files of each transect
FileList=dir('*ASC.TXT');

for n = 1:length(FileList)
    disp(['processing ' num2str(n)]);
    
    adcp(n) = adcp_ascii(FileList(n).name);
end

for jj=1:length(adcp)
    
        adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
        adcp(jj).hour,adcp(jj).minute,adcp(jj).second);
end

cd C:\GLOVER\output\myanmar\mmi_discharge
save('br_07mar2018','adcp')

%%
clear all,close all,clc
cd C:\GLOVER\output\myanmar\mmi_discharge

load('br_08mar2018')
int = adcp(1);
load('br_07mar2018')

int = [int;adcp(8);adcp(9)];
clear adcp
adcp = int;

for jj=1:3
    adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
        adcp(jj).hour,adcp(jj).minute,adcp(jj).second);
    t(jj)=adcp(jj).time(1);
end
[t,idx]=sort(t,'ascend');
adcp=adcp(idx);

save('br_int_march2018','adcp')

%% 
clear all%,close all,clc

cd C:\GLOVER\output\myanmar
load('mmi_discharge\br_int_march2018.mat')
load('longterminst\BogaleRiverInstruments.mat')
br=BogaleRiver;

for jj=1:3
adcp(jj).east_mean = nanmean(adcp(jj).east);
adcp(jj).north_mean = nanmean(adcp(jj).north);

adcp(jj).heading(adcp(jj).heading<0)=NaN;
adcp(jj).dir(adcp(jj).dir<0)=NaN;
adcp(jj).spd(adcp(jj).spd>200)=NaN;

adcp(jj).spd_mean = nanmean(adcp(jj).spd);
adcp(jj).spd_mean = movmean(adcp(jj).spd_mean,7);
adcp(jj).dir_mean = nanmean(adcp(jj).dir);
adcp(jj).dir_mean = movmean(adcp(jj).dir_mean,7);


figure(1);
subplot(211)
yyaxis left
plot(adcp(jj).time,adcp(jj).spd_mean,'k-'),hold on
yyaxis right
plot(br.datenum,br.ut_FredaDepth,'b-')
xlim([adcp(1).time(1) adcp(3).time(end)])
datetick('x','HH:MM','keeplimits')
% title(datestr(floor(adcp(jj).time(1))))
hold on
subplot(212)
yyaxis left
plot(adcp(jj).time,adcp(jj).dir_mean,'k-'),hold on
plot(adcp(jj).time,adcp(jj).heading,'r-'),hold on
yyaxis right
plot(br.datenum,br.ut_FredaDepth,'b-')
xlim([adcp(1).time(1) adcp(3).time(end)])
datetick('x','HH:MM','keeplimits')
hold on

% figure(2);
% scatter(adcp(jj).dir_mean,adcp(jj).heading,[],adcp(jj).time),hold on
% colorbar,refline(1,0),refline(-1,350)

adcp(jj).flowdiff = abs(adcp(jj).dir_mean - adcp(jj).heading);
figure(2);
scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).flowdiff,'.')
colorbar,caxis([0 180]),hold on,title('water spd')
colormap(cmocean('phase')),title('water dir - head')

figure;%(2)
scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).spd_mean)
colorbar,caxis([0 50]),hold on,title('water spd')
figure;%(3)
subplot(121),scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).dir_mean)
colorbar,caxis([0 360]),hold on,title('water dir')
subplot(122),scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).heading)
colorbar,caxis([0 360]), hold on
colormap(cmocean('phase')),title('head')
end


%%
clear all,close all,clc

cd C:\GLOVER\output\myanmar\mmi_discharge
load('br_int_march2018')

for jj=1:3
adcp(jj).dir(adcp(jj).dir<0)=NaN;
adcp(jj).spd(adcp(jj).spd>200)=NaN;

figure(10);
quiver(adcp(jj).lon,adcp(jj).lat,adcp(jj).east(2,:),adcp(jj).north(2,:))
hold on,

figure(20);
scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).spd(2,:))
hold on,colorbar,caxis([0 50]),title('water speed')

figure(30);
scatter(adcp(jj).lon,adcp(jj).lat,[],adcp(jj).depth(1,:)),colorbar
hold on,colorbar,title('depth')
end

%% add in CTD profiles of turbidity and salinity:

%% pull out ctd profiles from Int channel survey March 2018
clear all,close all,clc
cd C:\GLOVER\output\myanmar

load('AyeMar18_CTD_all.mat')

for jj=1:length(CombinedProfiles)
    
    t(jj)=floor(CombinedProfiles(jj).time(1));
end
[t,idx]=sort(t,'ascend');
CombinedProfiles=CombinedProfiles(idx);
ctd=CombinedProfiles(t==datenum('8-March-2018'));

save('C:\GLOVER\output\myanmar\mmi_discharge\ctdprof_8Mar18.mat','ctd')

%% pull out ctd profiles from Int channel survey Sep 2017
clear all,close all,clc
cd C:\GLOVER\output\myanmar

load('AyeSept17_CTD_all.mat')

for jj=1:length(CombinedProfiles)
    
    t(jj)=floor(CombinedProfiles(jj).time(1));
end
[t,idx]=sort(t,'ascend');
CombinedProfiles=CombinedProfiles(idx);
ctd=CombinedProfiles(t==datenum('13-Sept-2017'));

save('C:\GLOVER\output\myanmar\mmi_discharge\ctdprof_13Sep17.mat','ctd')

%% add dist along chnl to ctd profile struct
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('mmi_discharge\ctdprof_13Sep17.mat')
load('mmi_discharge\channel_trace.mat')

%figure out where ctd profs are along channel

for jj=1:length(ctd)
    [x,y,~] = deg2utm(ctd(jj).lat, ctd(jj).long);
    dst = sqrt(abs(ch.x-x).^2 + abs(ch.y-y).^2);
    [~,idx]=min(dst);
    ctd(jj).dist=ch.dist_in(idx);
end