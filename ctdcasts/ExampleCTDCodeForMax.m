clear all,close all,clc

load('AyeMar18_CTD_all.mat')

%% Pull out all the casts that happened on a specific date:

% Get the day of month for each cast:
for jj=1:length(CombinedProfiles)
    daychk(jj,:)=datevec(CombinedProfiles(jj).time(1));
end
daychk=daychk(:,3);

% sort the casts in order by day of month:
[daychk,idx]=sort(daychk);
CombinedProfiles=CombinedProfiles(idx);
% remove repeat profiles and mangrove forest profs
CombinedProfiles([9,12,17,95:101])=[];
daychk([9,12,17,95:101])=[];

% make new structure that is jsut the profiles from 1 day:
dd = CombinedProfiles(daychk==14);

%% Create a single matrix with depth on y, and distance along river on x
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
    md(jj)=max(dd(jj).SSCCal);
    max(dd(jj).Salinity)
end
% define new distances along the idealized transect 0=RK58
dist =  round(sqrt(prj(:,1).^2 + prj(:,2).^2)/100); % 


%% Make empty matrix for ssc/sal at 0.05cm depth interval, for all points
% along the river distance
L=370; %length upriver in m/100
depthvec=[0:0.05:18]';
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
% fill surface (upper 3cm) values using nearest
sal_vec(1:60,:)=fillmissing(sal_vec(1:60,:),'nearest');
ssc_vec(1:60,:)=fillmissing(ssc_vec(1:60,:),'nearest');

%%
% fill laterally using linear interp
maxdepth=fillmissing(maxdepth,'linear');
sal_vec=fillmissing(sal_vec,'linear',2);%linear or nearest
ssc_vec=fillmissing(ssc_vec,'linear',2);

% remove vals below max depth
for jj=1:length(maxdepth)
    idx=find(depthvec>maxdepth(jj));
    sal_vec(idx,jj)=NaN;
    ssc_vec(idx,jj)=NaN;
end
% remove vals seaward of first cast
%sal_vec(:,1:min(dist)-10)=NaN;
%ssc_vec(:,1:min(dist)-10)=NaN;
sal_vec(sal_vec<0)=NaN;
ssc_vec(ssc_vec<0 | ssc_vec>870)=NaN;

%% Plot the figure
% define colormap:
C1=linspace(0.6,0.95,400)';
C1(:,2)=linspace(0.1,0.95,400)';
C1(:,3)=linspace(0.1,1,400)';
C1=flipud(C1);
% plot color as SSC with lines of Sal overlying
levels=0:2:15;

subplot(313)
pcolor(1:L,depthvec,ssc_vec),shading interp,axis ij,hold on
colorbar,colormap(C1),caxis([0 1000])
% contour(1:L,depthvec,sal_vec,levels,'showtext','on','Color','k'),hold on
plot(dist,zeros(1,length(dist)),'kv')
ylim([0 21]),axis ij

