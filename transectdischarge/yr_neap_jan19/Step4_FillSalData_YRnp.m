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
% units of mg/L

clear all,close all,clc
load('AyeJan19_CTD_all.mat')
fid=fopen('SSC_ForProfiles.csv');
T=textscan(fid,'%f %f','Delimiter',',');fclose(fid);

sal=profiles(1:12); clear profiles

for jj=1:length(sal)
    sal(jj).timeidx=nanmean(sal(jj).time);%time cast starts, for matching adcp transect
    sal(jj).Salinity(end-4:end)=sal(jj).Salinity(end-4);
end

load('YR_Jan19_Neap.mat')

allidx_adcp=[2,4,6,8,10,11,12,13,15,17,19,22];
for jj=1:length(sal)
    % get the adcp transect
    sal(jj).transect=allidx_adcp(jj);
end

%% find idealized transect end points (from adcp processing:
% BR_neap_Residuals)
[right.x,right.y]=deg2utm(16.758914, 96.168642);
[left.x,left.y]=deg2utm(16.766158, 96.170355);
for jj=1:length(sal)   
    % convert the coordinates of your ADCP records to UTM
    [sal(jj).x, sal(jj).y,~] = deg2utm(sal(jj).lat, sal(jj).long);
    % project the ADCP data onto the idealized transect
    salproj = proj([left.x - sal(jj).x, left.y - sal(jj).y],...
        [left.x- right.x, left.y - right.y]);
    % define new distances along the idealized transect
    sal(jj).dist =  round(sqrt(salproj(1).^2 + salproj(2).^2));
end
% now you have indexed CTD profiles, put those into useful grid matching
% the adcp style:
% 0:1000 wide (1 m bin width)
[~,cols]=size(adcp(1).alongComplete);
% 0:0.01:20m deep (change to exceed max depth of adcp data)
interpdepth=0:0.05:20;

%%
for jj=1:length(sal)
    sal(jj).interpalong=NaN([length(interpdepth),cols]);
    % add the sal profile to the correct location along the transect and
    % interpolate it onto a consistent depth matrix
    sal(jj).interpalong(:,sal(jj).dist)=interp1(sal(jj).Depth,sal(jj).Salinity,interpdepth);
    
    % fill the surface and bed with avg of nearest value
    fillsurface=find(interpdepth<min(sal(jj).Depth));
    sal(jj).interpalong(fillsurface,sal(jj).dist)=nanmean(...
        sal(jj).interpalong(fillsurface(end):fillsurface(end)+3,sal(jj).dist));
    fillbed=find(interpdepth>max(sal(jj).Depth));
    sal(jj).interpalong(fillbed,sal(jj).dist)=nanmean(...
        sal(jj).interpalong(fillbed(1)-10:fillbed(1),sal(jj).dist));
    sal(jj).adcpbeddepth=adcp(sal(jj).transect).interpdepths;
    

clear fill*
end
clear *adcp iC iA T

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


save('YR_Jan19_Neap_Sal','sal')