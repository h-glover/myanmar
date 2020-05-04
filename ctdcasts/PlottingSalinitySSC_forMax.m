clear all,close all,clc
% define colormap:
C1=linspace(0.6,0.95,400)';
C1(:,2)=linspace(0.1,0.95,400)';
C1(:,3)=linspace(0.1,1,400)';
C1=flipud(C1);

load('AyeMar18_CTD_all.mat')
for jj=1:length(CombinedProfiles)
    daychk(jj,:)=datevec(CombinedProfiles(jj).time(1));
end
daychk=daychk(:,3);
[daychk,idx]=sort(daychk);
CombinedProfiles=CombinedProfiles(idx);
CombinedProfiles([9,12,17,95:101])=[];
daychk([9,12,17,95:101])=[];% remove repeat profiles and mangrove forest profs


% Yangon
dd = CombinedProfiles(daychk==1 | daychk==2);

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
dist =  sqrt(prj(:,1).^2 + prj(:,2).^2)/1000;
[dist,idx]=sort(dist,'ascend');
dd=dd(idx);

% Make empty matrix for ssc/sal at 0.05cm depth interval, for all points
% along the river distance
L=length(dist); %length upriver in m/100
depthvec=[0:0.05:27]';
sal_vec=NaN(length(depthvec),L);
ssc_vec=NaN(length(depthvec),L);

for jj=1:length(dd)
    dd(jj).Depth=fillmissing(dd(jj).Depth,'linear');
    ssc_vec(:,jj)=interp1(dd(jj).Depth,dd(jj).SSCCal,depthvec);
    sal_vec(:,jj)=interp1(dd(jj).Depth,dd(jj).Salinity,depthvec);
end
% fill surface (upper 3m) values using nearest
ssc_vec(1:60,:)=fillmissing(ssc_vec(1:60,:),'nearest');
sal_vec(1:60,:)=fillmissing(sal_vec(1:60,:),'nearest');

% plot color as SSC with lines of Salinity overlying at 2 unit intervals
levels=0:2:15;
figure;
pcolor(dist,depthvec,ssc_vec),shading interp,axis ij,hold on
colorbar,colormap(C1),caxis([0 4000]) %change the color axis for Bogale and Pathein
contour(dist,depthvec,sal_vec,levels,'showtext','on','Color','k')
hold on,axis ij
ylim([0 21]),xlim([0 55])

