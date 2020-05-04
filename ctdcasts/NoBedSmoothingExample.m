%%
clear all,%close all,clc
load('AyeSept17_CTD_all.mat')
% define colormap:
C1=linspace(0.6,0.95,400)';
C1(:,2)=linspace(0.1,0.95,400)';
C1(:,3)=linspace(0.1,1,400)';
C1=flipud(C1);

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
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2));
[dist,srt]= sort(dist,'ascend');
dd=dd(srt);

% Make empty matrix for ssc/sal at 0.05cm depth interval, for all points
% along the river distance
L=length(dd); %length upriver in m/100
depthvec=[0:0.05:27]';
sal_vec=NaN(length(depthvec),L);
ssc_vec=NaN(length(depthvec),L);
% find river depth so vals can be removed below this depth
maxdepth=NaN(1,L);

for jj=1:L
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    ssc_vec(:,jj)=interp1(dd(jj).Depth,dd(jj).SSCCal,depthvec);
    maxdepth(1,jj)=max(dd(jj).Depth);
end
% fill surface (upper 3m) values using nearest
ssc_vec(1:60,:)=fillmissing(ssc_vec(1:60,:),'nearest');
% fill laterally using linear interp
maxdepth=fillmissing(maxdepth,'pchip');
% maxdepth(maxdepth>21)=21;

% remove vals below max depth
for jj=1:length(maxdepth)
    idx=find(depthvec>maxdepth(jj));
    ssc_vec(idx,jj)=NaN;
end
ssc_vec(ssc_vec<0)=NaN;

% plot color as SSC with lines of Sal overlying
levels=0:1:15;
figure;
subplot(311)
pcolor(dist/1000,depthvec,ssc_vec),shading interp,axis ij,hold on
colorbar,colormap(C1),caxis([0 700])
plot(dist/1000,zeros(1,length(dist)),'kv'),axis ij
ylim([0 21])
xlim([5 55])

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
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2));
[dist,srt]= sort(dist,'ascend');
dd=dd(srt);

% Make empty matrix for ssc/sal at 0.05cm depth interval, for all points
% along the river distance
L=length(dd); %length upriver in m/100
depthvec=[0:0.05:27]';
sal_vec=NaN(length(depthvec),L);
ssc_vec=NaN(length(depthvec),L);
% find river depth so vals can be removed below this depth
maxdepth=NaN(1,L);

for jj=1:L
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    ssc_vec(:,jj)=interp1(dd(jj).Depth,dd(jj).SSCCal,depthvec);
    sal_vec(:,jj)=interp1(dd(jj).Depth,dd(jj).Salinity,depthvec);
    maxdepth(1,jj)=max(dd(jj).Depth);
end
% fill surface (upper 3m) values using nearest
ssc_vec(1:60,:)=fillmissing(ssc_vec(1:60,:),'nearest');
% fill laterally using linear interp
maxdepth=fillmissing(maxdepth,'pchip');
% maxdepth(maxdepth>21)=21;

% remove vals below max depth
for jj=1:length(maxdepth)
    idx=find(depthvec>maxdepth(jj));
    ssc_vec(idx,jj)=NaN;
end
ssc_vec(ssc_vec<0)=NaN;

levels=0:2:10;
subplot(312)
pcolor(dist/1000,depthvec,ssc_vec),shading interp,axis ij,hold on
colorbar,colormap(C1),caxis([0 700])
contour(dist/1000,depthvec,sal_vec,levels,'showtext','on','Color','k')
axis ij,hold on
plot(dist/1000,zeros(1,length(dist)),'kv')
ylim([0 21])
xlim([5 55])

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
% define new distances along the idealized transect 0=RK58
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2))+58000;
[dist,srt]= sort(dist,'ascend');
dd=dd(srt);

% Make empty matrix for ssc/sal at 0.05cm depth interval, for all points
% along the river distance
L=length(dd); %length upriver in m/100
depthvec=[0:0.05:27]';
sal_vec=NaN(length(depthvec),L);
ssc_vec=NaN(length(depthvec),L);
% find river depth so vals can be removed below this depth
maxdepth=NaN(1,L);

for jj=1:L
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    ssc_vec(:,jj)=interp1(dd(jj).Depth,dd(jj).SSCCal,depthvec);
    maxdepth(1,jj)=max(dd(jj).Depth);
end
% fill surface (upper 3m) values using nearest
ssc_vec(1:60,:)=fillmissing(ssc_vec(1:60,:),'nearest');
% fill laterally using linear interp
maxdepth=fillmissing(maxdepth,'pchip');
% maxdepth(maxdepth>21)=21;

% remove vals below max depth
for jj=1:length(maxdepth)
    idx=find(depthvec>maxdepth(jj));
    ssc_vec(idx,jj)=NaN;
end
ssc_vec(ssc_vec<0)=NaN;


subplot(313)
pcolor(dist/1000,depthvec,ssc_vec),shading interp,axis ij,hold on
colorbar,colormap(C1),caxis([0 700])
plot(dist/1000,zeros(1,length(dist)),'kv'),axis ij
ylim([0 21])
% xlim([0 ])