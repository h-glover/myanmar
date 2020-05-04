%
clear all

load('BR_Sept17_Spring.mat')
load('AyeSept17_CTD_all.mat')
%
ctd=CombinedProfiles; clear CombinedProfiles
for jj=1:length(ctd)
    ctd(jj).profileidx=jj; %label 1:length(ctd) for later use
    ctd(jj).timeidx=ctd(jj).time(1);%time cast starts, for matching adcp transect
end
% concatenate the profile indexs and times
allidx_ctd=horzcat(ctd.profileidx);
alltime_ctd=vertcat(ctd.timeidx);
alltime_ctd=datevec(alltime_ctd);alltime_ctd(:,6)=0;
alltime_ctd=datenum(alltime_ctd);

for jj=1:length(adcp)
    %get the matlab time for each transect and assign a transect number to
    %each time point
    adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
        adcp(jj).hour+6,adcp(jj).minute,zeros(1,length(adcp(jj).minute)));
    adcp(jj).transidx=jj*ones(1,length(adcp(jj).time));
end
% concatenate transect index into 1 vector 
allidx_adcp=horzcat(adcp.transidx);
alltime_adcp=horzcat(adcp.time);
% compare adcp time and ctd time (minutes removed)
[~,iA,iC]=intersect(alltime_adcp,alltime_ctd);
% pull out all ctd profiles from this day
ctd_trans=allidx_ctd(iC);
allidx_adcp=allidx_adcp(iA);
ctd=ctd(iC);
% find idealized transect end points (from adcp processing:
% BR_neap_Residuals)
[right.x,right.y]=deg2utm(16.099239, 95.320838);
[left.x,left.y]=deg2utm(16.100713, 95.330033);
for jj=1:length(ctd)
    % get the adcp transect
    ctd(jj).transect=allidx_adcp(jj);
    % convert the coordinates of your ADCP records to UTM
    [ctd(jj).x, ctd(jj).y,~] = deg2utm(ctd(jj).lat, ctd(jj).long);
    % project the ADCP data onto the idealized transect
    ctdproj = proj([left.x - ctd(jj).x, left.y - ctd(jj).y],...
        [left.x- right.x, left.y - right.y]);
    % define new distances along the idealized transect
    ctd(jj).dist =  round(sqrt(ctdproj(1).^2 + ctdproj(2).^2));
end
% 
for jj=1:length(ctd)
    if ctd(jj).dist<200
        stnidx(jj)=1;
    elseif ctd(jj).dist>700
        stnidx(jj)=3;
    else
        stnidx(jj)=2;
    end
end

clear all* i* ctd_trans left right ctdproj

depthvec=linspace(0,16.2,100);
salvec=NaN([length(depthvec),length(ctd)]);
sscvec=NaN([length(depthvec),length(ctd)]);
for jj=1:length(ctd)
    idx=find(~isnan(ctd(jj).Depth));
    salvec(:,jj)=interp1(ctd(jj).Depth(idx),ctd(jj).Salinity(idx),depthvec);
    sscvec(:,jj)=interp1(ctd(jj).Depth(idx),ctd(jj).SSCCal(idx),depthvec);
    timevec(jj)=ctd(jj).time(1);
end
%save('BR_neap_Salinity','salvec','timevec','depthvec','stnidx')




%%
close all
load('BR_Sept17_Spring_SedFlux.mat')

depths=0:0.5:30.5; %change based on length of adcp alongComplete

C1=flipud([ones([30,1]),linspace(1,0,30)',linspace(1,0,30)']);
C2=[linspace(1,0,30)',linspace(1,0,30)',ones([30,1])];
C=[C1(1:end-3,:);C2(4:end,:)];
for jj=[9,10,32,33]
    F=figure;
    F.Position=[360.3333 265.6667 824.6667 352.0000];
    p=pcolor(0:1000,depths,adcp(jj).alongComplete);shading interp,axis ij
    p.EdgeColor='interp';
    c=colorbar;
    colormap(C)
    caxis([-1 1])
    title(['Transect ',num2str(jj),' ',datestr(adcp(jj).time(1),'mm/dd/yy HH:MM')])
    ylabel('')
    c.Label.String='m/s';
    hold on
    ylim([0 16]),xlim([50 1000])
end

%% SSC profiles to go with above transect examples
close all
transidx=vertcat(ctd.transect);
%[9,10,32,33]
figure;
subplot(121)
profidx=find(transidx==13);
for jj=1:length(profidx)
    plot(ctd(profidx(jj)).SSCCal(4:end),ctd(profidx(jj)).Depth(4:end),'k')
    hold on
end
axis ij,xlim([0 600]),ylim([0 16])
xlabel('SSC (mg/L)'),ylabel('Depth (m)')

subplot(122)
profidx=find(transidx==30 | transidx==31);
for jj=1:length(profidx)
    plot(ctd(profidx(jj)).SSCCal(4:end),ctd(profidx(jj)).Depth(4:end),'k')
    hold on
end
axis ij,xlim([0 600]),ylim([0 16])
xlabel('SSC (mg/L)'),ylabel('Depth (m)')