clear all
cd C:\GLOVER\output\myanmar\
%%
clear all,%close all,clc
load('AyeSept17_CTD_all.mat')

C1=cmocean('turbid');C1=C1(1:200,:);

for jj=1:length(CombinedProfiles)
    daychk(jj,:)=datevec(CombinedProfiles(jj).time(1));
end
daychk=daychk(:,3);
[daychk,idx]=sort(daychk);
CombinedProfiles=CombinedProfiles(idx);

dd = CombinedProfiles(daychk==17);

% project the ctd locations to the thalweg to get distance along estuary
lat=vertcat(dd.lat);
lon=vertcat(dd.long);
[loc.x,loc.y]=deg2utm(lat,lon);

% define 3 transects along estuary for distance along river calculation
[mouth.x,mouth.y]=deg2utm(16.447348, 96.352243);
[split.x,split.y]=deg2utm( 16.756927, 96.199459);
[turn.x,turn.y]=deg2utm( 16.775776, 96.119558);
[top.x,top.y]=deg2utm(16.880264,96.087471);
% project the casts onto the nearest transect along the river
for jj=1:length(dd)
    if loc.x(jj)>=split.x %lower river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);
    elseif loc.x(jj)>=turn.x && loc.x(jj)<split.x %middle river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);    
    elseif loc.x(jj)<turn.x %upper river
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- split.x, mouth.y- split.y]);    
    end
end
% define new distances along the idealized transect 0=mouth
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2)/100);


% Make empty matrix for ssc/sal at 0.05cm depth interval, for all points
% along the river distance
L=570; %length upriver in m/100
depthvec=[0:0.05:27]';
sal_vec=NaN(length(depthvec),L);
ssc_vec=NaN(length(depthvec),L);
% find river depth so vals can be removed below this depth
maxdepth=NaN(1,L);
cidx=(1/5)*ones(5,1);
for jj=1:length(dd)
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    ssc_vec(:,dist(jj))=interp1(dd(jj).Depth,dd(jj).SSCCal,depthvec);
    ssc_vec(:,dist(jj))=conv(ssc_vec(:,dist(jj)),cidx,'same');
    maxdepth(1,dist(jj))=max(dd(jj).Depth);
end
% fill surface (upper 3m) values using nearest
ssc_vec(1:60,:)=fillmissing(ssc_vec(1:60,:),'nearest');
% fill laterally using linear interp
maxdepth=fillmissing(maxdepth,'nearest');
% maxdepth(maxdepth>21)=21;
% figure;plot(maxdepth)
ssc_vec=fillmissing(ssc_vec,'linear',2);
K = (1/15)*ones(5);
%ssc_vec= conv2(ssc_vec,K,'same');
% remove vals below max depth
for jj=1:length(maxdepth)
    idx=find(depthvec>maxdepth(jj));
    ssc_vec(idx,jj)=NaN;
end
% remove vals seaward of first cast
ssc_vec(:,1:min(dist)-10)=NaN;
ssc_vec(:,max(dist)+10:end)=NaN;
ssc_vec(ssc_vec<0)=NaN;

% plot color as SSC with lines of Sal overlying
levels=0:1:15;
figure;
subplot(311)
pcolor(1:L,depthvec,ssc_vec),shading interp,axis ij,hold on
colorbar,colormap(C1),caxis([0 700])
% plot(dist,zeros(1,length(dist)),'kv'),axis ij
ylim([0 21])
xlim([50 550])

%% Bogale
clearvars -except CombinedProfiles daychk C1
dd = CombinedProfiles(daychk==12);

% project the ctd locations to the thalweg to get distance along estuary
lat=vertcat(dd.lat);
lon=vertcat(dd.long);
[loc.x,loc.y]=deg2utm(lat,lon);

[mouth.x,mouth.y]=deg2utm(15.755580, 95.229443);
[top.x,top.y]=deg2utm( 16.118759, 95.324648);
% project the casts onto the nearest transect along the river
for jj=1:length(dd)
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- top.x, mouth.y- top.y]);
    md(jj)=max(dd(jj).Depth);
        ms(jj)=max(dd(jj).SSCCal);
end
% define new distances along the idealized transect 0=RK58
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2)/100); % 

% Make empty matrix for ssc/sal at 0.05cm depth interval, for all points
% along the river distance
L=400; %length upriver in m/100
depthvec=[0:0.05:13]';
sal_vec=NaN(length(depthvec),L);
ssc_vec=NaN(length(depthvec),L);
% find river depth so vals can be removed below this depth
maxdepth=NaN(1,L);
for jj=1:length(dd)
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    
    sal_vec(:,dist(jj))=interp1(dd(jj).Depth,dd(jj).Salinity,depthvec);
    ssc_vec(:,dist(jj))=interp1(dd(jj).Depth,dd(jj).SSCCal,depthvec);
    maxdepth(1,dist(jj))=max(dd(jj).Depth);
end
% fill surface (upper 3m) values using nearest
sal_vec(1:60,:)=fillmissing(sal_vec(1:60,:),'nearest');
ssc_vec(1:60,:)=fillmissing(ssc_vec(1:60,:),'nearest');
% fill laterally using linear interp
maxdepth=fillmissing(maxdepth,'nearest');
% figure;plot(maxdepth)
sal_vec=fillmissing(sal_vec,'linear',2);
ssc_vec=fillmissing(ssc_vec,'linear',2);
% remove vals below max depth
for jj=1:length(maxdepth)
    idx=find(depthvec>maxdepth(jj));
    sal_vec(idx,jj)=NaN;
    ssc_vec(idx,jj)=NaN;
end
% remove vals seaward of first cast
sal_vec(:,1:min(dist)-10)=NaN;
ssc_vec(:,1:min(dist)-10)=NaN;
sal_vec(sal_vec<0)=NaN;
ssc_vec(ssc_vec<0)=NaN;

% plot color as SSC with lines of Sal overlying
levels=0:2:10;

subplot(312)
pcolor(1:L,depthvec,ssc_vec),shading interp,axis ij,hold on
colorbar,caxis([0 700])
% contour(1:L,depthvec,sal_vec,levels,'showtext','on','Color','k')
% axis ij,hold on
% plot(dist,zeros(1,length(dist)),'kv')
ylim([0 21])
xlim([50 550])
%% Pathein
clearvars -except CombinedProfiles daychk C1
dd = CombinedProfiles(daychk==20);

% project the ctd locations to the thalweg to get distance along estuary
lat=vertcat(dd.lat);
lon=vertcat(dd.long);
[loc.x,loc.y]=deg2utm(lat,lon);

[mouth.x,mouth.y]=deg2utm(16.354683, 94.676328);
[top.x,top.y]=deg2utm(16.778223,94.728839);
% project the casts onto the nearest transect along the river
for jj=1:length(dd)
        prj(jj,:) = proj([mouth.x - loc.x(jj),...
            mouth.y - loc.y(jj)], [mouth.x- top.x, mouth.y- top.y]);
    md(jj)=max(dd(jj).Depth);ms(jj)=max(dd(jj).SSCCal);
end
% define new distances along the idealized transect 0=RK58
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2)/100); % 


% Make empty matrix for ssc/sal at 0.05cm depth interval, for all points
% along the river distance
L=370; %length upriver in m/100
depthvec=[0:0.05:25]';
ssc_vec=NaN(length(depthvec),L);
% find river depth so vals can be removed below this depth
maxdepth=NaN(1,L);
for jj=1:length(dd)
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    
    ssc_vec(:,dist(jj))=interp1(dd(jj).Depth,dd(jj).SSCCal,depthvec);
    maxdepth(1,dist(jj))=max(dd(jj).Depth);
end
% fill surface (upper 3m) values using nearest
ssc_vec(1:60,:)=fillmissing(ssc_vec(1:60,:),'nearest');
% fill laterally using linear interp
maxdepth=fillmissing(maxdepth,'nearest');
% figure;plot(maxdepth)
ssc_vec=fillmissing(ssc_vec,'linear',2);
% remove vals below max depth
for jj=1:length(maxdepth)
    idx=find(depthvec>maxdepth(jj));
    ssc_vec(idx,jj)=NaN;
end
% remove vals seaward of first cast
ssc_vec(:,[1:min(dist)-10,max(dist)+10:end])=NaN;
ssc_vec(ssc_vec<0)=NaN;

% plot color as SSC with lines of Sal overlying
levels=0:2:15;

subplot(313)
pcolor(1:L,depthvec,ssc_vec),shading interp,axis ij,hold on
colorbar,caxis([0 700])
% plot(dist,zeros(1,length(dist)),'kv'),axis ij
ylim([0 21])
xlim([0 370])