% Use this code to convert individual CTD casts from the compiled list into
% a series of interpolated transects associated with an adcp transect. This
% code first compiles the casts by transect, then converts to sigma
% coordinates, then interpolates.
% 1. Identify casts in the same transect
% 2. put casts onto an individual transect matrix: X rows deep based on
% deepest cast and y columns wide based on standard river width.
% 3. Project casts onto transect based on GPS pt to get along transect location
% 4. Pad this matrix with values to seabed
% 5. Convert to sigma coordinates
% 6. Interpolate across transects

clear all,close all,clc
load('AyeMar18_CTD_all.mat')

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


load('BR_Mar18_Spring.mat')
% % adcp time and index will be reused in the SedFlux calculation
% % only run this section once, the first time you process
% for jj=1:length(adcp)
%     %get the matlab time for each transect and assign a transect number to
%     %each time point
%     adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
%         adcp(jj).hour,adcp(jj).minute,zeros(1,length(adcp(jj).minute)));
%     adcp(jj).transidx=jj*ones(1,length(adcp(jj).time));
% end
% save('BR_Mar18_Spring','adcp') 

% concatenate transect index into 1 vector 
allidx_adcp=horzcat(adcp.transidx);
alltime_adcp=horzcat(adcp.time);
% compare adcp time and ctd time (minutes removed)
[~,iA,iC]=intersect(alltime_adcp,alltime_ctd);
% pull out all ctd profiles from this day
ctd_trans=allidx_ctd(iC);
allidx_adcp=allidx_adcp(iA);
ctd=ctd(iC);
for jj=1:length(ctd)
    % get the adcp transect
    ctd(jj).transect=allidx_adcp(jj);
    L(jj)=length(ctd(jj).time);
end
rows=max(L);

% find distance along transect
% find idealized transect end points (from adcp processing:
% BR_neap_Residuals)
[right.x,right.y]=deg2utm(16.099239, 95.320838);
[left.x,left.y]=deg2utm(16.100713, 95.330033);
for jj=1:length(ctd)   
    % convert the coordinates of your ADCP records to UTM
    [ctd(jj).x, ctd(jj).y,~] = deg2utm(ctd(jj).lat, ctd(jj).long);
    % project the ADCP data onto the idealized transect
    ctdproj = proj([left.x - ctd(jj).x, left.y - ctd(jj).y],...
        [left.x- right.x, left.y - right.y]);
    % define new distances along the idealized transect
    ctd(jj).dist =  round(sqrt(ctdproj(1).^2 + ctdproj(2).^2));
end


% now you have indexed CTD profiles, put those into useful grid matching
% the adcp style:
% 0:1000 wide (1 m bin width)
[~,cols]=size(adcp(1).alongComplete);
% 0:0.01:20m deep (change to exceed max depth of adcp data)
interpdepth=0:0.05:20;

% make sal structure for first transect CTD
% then for 2:length(ctd), either:
% add to the current matrix
% or make a new one if its from different transect than the previous
L=length(ctd(1).Depth);
ctd(jj).Depth=fillmissing(ctd(jj).Depth,'linear');
sal.transect=ctd(1).transect;
sal.lat=ctd(1).lat;
sal.long=ctd(1).long;
sal.dist=ctd(1).dist;
sal.Density_avg=ctd(jj).Density_avg;
sal.time=NaN([rows,3]);
sal.time(1:L,1)=ctd(1).time;
sal.Depth=NaN([rows,3]);
sal.Depth(1:L,1)=ctd(1).Depth;
sal.sal=NaN([rows,3]);
sal.sal(1:L,1)=ctd(1).Salinity;

% add the sal profile to the correct location along the transect and
% interpolate it onto a consistent depth matrix
sal.interpalong=NaN([length(interpdepth),cols]);
sal.interpalong(:,ctd(1).dist)=interp1(ctd(1).Depth,ctd(1).Salinity,interpdepth);

% fill the surface and bed with avg of nearest value
fillsurface=find(interpdepth<min(ctd(1).Depth));
sal.interpalong(fillsurface,ctd(1).dist)=nanmean(...
    sal.interpalong(fillsurface(end):fillsurface(end)+3,ctd(1).dist));
fillbed=find(interpdepth>max(ctd(1).Depth));
sal.interpalong(fillbed,ctd(1).dist)=nanmean(...
    sal.interpalong(fillbed(1)-3:fillbed(1),ctd(1).dist));

% add in the adcp bed depth for that transect (used for sigma conversion)
sal.adcpbeddepth=adcp(sal(1).transect).interpdepths;

%if next profile is in the same transect then add data in a new column in the same matrix
col=2; 
% counter is the sal() number (number of transects)
counter=1;
for jj=2:length(ctd)
    L=length(ctd(jj).Depth);
    ctd(jj).Depth=fillmissing(ctd(jj).Depth,'linear');
    if ctd(jj).transect==ctd(jj-1).transect %transect is same as prev cast
        sal(counter).lat(col)=ctd(jj).lat;
        sal(counter).long(col)=ctd(jj).long;
        sal(counter).dist(col)=ctd(jj).dist;
        sal(counter).Density_avg(col)=ctd(jj).Density_avg;
        sal(counter).time(1:L,col)=ctd(jj).time;
        sal(counter).Depth(1:L,col)=ctd(jj).Depth;
        sal(counter).sal(1:L,col)=ctd(jj).Salinity;
        col=col+1; %add new column for the next profile
    else
        counter=counter+1; %add new sal row to structure
        sal(counter).transect=ctd(jj).transect;
        sal(counter).lat(1)=ctd(jj).lat;
        sal(counter).long(1)=ctd(jj).long;
        sal(counter).dist(1)=ctd(jj).dist;
        sal(counter).Density_avg(1)=ctd(jj).Density_avg;
        sal(counter).maxdepth(1)=max(ctd(jj).Depth);
        sal(counter).time=NaN([rows,3]);
        sal(counter).time(1:L,1)=ctd(jj).time;
        sal(counter).Depth=NaN([rows,3]);
        sal(counter).Depth(1:L,1)=ctd(jj).Depth;
        sal(counter).sal=NaN([rows,3]);
        sal(counter).sal(1:L,1)=ctd(jj).Salinity;
        sal(counter).interpalong=NaN([length(interpdepth),cols]); %rows=0.0001m bins, cols=same as ADCP, 1m bins
        % get transect depths from adcp transect:
        sal(counter).adcpbeddepth=adcp(sal(counter).transect).interpdepths;
        col=2;
    end
    % add the sal profile to the correct location along the transect and
    % interpolate it onto a consistent depth matrix
    sal(counter).interpalong(:,ctd(jj).dist)=interp1(ctd(jj).Depth,ctd(jj).Salinity,interpdepth);
    
    % fill the surface and bed with avg of nearest value
    fillsurface=find(interpdepth<min(ctd(jj).Depth));
    sal(counter).interpalong(fillsurface,ctd(jj).dist)=nanmean(...
        sal(counter).interpalong(fillsurface(end):fillsurface(end)+3,ctd(jj).dist));
    fillbed=find(interpdepth>max(ctd(jj).Depth));
    sal(counter).interpalong(fillbed,ctd(jj).dist)=nanmean(...
        sal(counter).interpalong(fillbed(1)-3:fillbed(1),ctd(jj).dist));
    
    clear fill*
end
clear *adcp c* iA iC

%% If there's only 1 cast, use the previous cast data to fill out the
%transect
for jj=2:length(sal)
    if length(sal(jj).dist)==1 && length(sal(jj-1).dist)==3
        sal(jj).dist(:,2:3)=sal(jj-1).dist(:,2:3);
        sal(jj).lat(:,2:3)=sal(jj-1).lat(:,2:3);
        sal(jj).long(:,2:3)=sal(jj-1).long(:,2:3);
        sal(jj).time(:,2:3)=sal(jj-1).time(:,2:3);
        sal(jj).Depth(:,2:3)=sal(jj-1).Depth(:,2:3);
        sal(jj).sal(:,2:3)=sal(jj-1).sal(:,2:3);
        sal(jj).interpalong(:,sal(jj).dist(1)+1:end)=sal(jj-1).interpalong(:,(sal(jj).dist(1)+1):end);
        
    elseif length(sal(jj).dist)==1 && length(sal(jj-1).dist)==2
        sal(jj).dist(:,2)=sal(jj-1).dist(:,2);
        sal(jj).lat(:,2)=sal(jj-1).lat(:,2);
        sal(jj).long(:,2)=sal(jj-1).long(:,2);
        sal(jj).time(:,2)=sal(jj-1).time(:,2);
        sal(jj).Depth(:,2)=sal(jj-1).Depth(:,2);
        sal(jj).sal(:,2)=sal(jj-1).sal(:,2);
        sal(jj).interpalong(:,sal(jj).dist(2))=sal(jj-1).interpalong(:,sal(jj).dist(2));
    end
    
end

%% Put the ctd profiles into sigma coordinates in a new matrix:
% 100pt vector between 0 and 1 for the sigma depths
% calc the sigma depths for each column (divide depth by max depth)
% interpolate onto the 100pt matrix
sigmadepth=linspace(0,1,50);

for jj=1:length(sal)
    % fill to bed with nearest value
    sal(jj).interpalong=fillmissing(sal(jj).interpalong,'nearest',1);
    
    % divide depth by max depth:
    sal(jj).sigmaDepth=interpdepth'./sal(jj).adcpbeddepth;
    % make a nan matrix that will be the sigma interpolation matrix
    sal(jj).sigmaAlongComplete=NaN([length(sigmadepth),1001]);
    
    for nn=sal(jj).dist
        sal(jj).sigmaAlongComplete(:,nn)=interp1(sal(jj).sigmaDepth(:,nn)...
            ,sal(jj).interpalong(:,nn),sigmadepth);
    end
end


%% next interpolate sal profiles across transects:

for jj=1:length(sal)
    % very simple QC/data clean up: sal should be positive
    sal(jj).sigmaAlongComplete(sal(jj).sigmaAlongComplete<=0)=NaN;
    
    % fill linearly across river
    sal(jj).sigmaAlongComplete=fillmissing(sal(jj).sigmaAlongComplete,'nearest',2);
    
    figure;
    pcolor(0:1000,sigmadepth,sal(jj).sigmaAlongComplete)
    shading flat,axis ij,colorbar
end



save('BR_Mar18_Spring_Sal','sal')