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

ssc=profiles(1:12); clear profiles

for jj=1:length(ssc)
    ssc(jj).timeidx=nanmean(ssc(jj).time);%time cast starts, for matching adcp transect
    
    % add in the near-bed ssc samples
    if T{2}(jj)==0
        ssc(jj).NBssc=ssc(jj-1).NBssc;
        ssc(jj).Surfssc=ssc(jj-1).Surfssc;
    else
        ssc(jj).NBssc=T{2}(jj).*1000;
        ssc(jj).Surfssc=T{1}(jj).*1000;
    end
    bed=max(ssc(jj).Depth);
    NBssc_prof=nanmean(ssc(jj).SSC(ssc(jj).Depth>(bed-3)));
    if NBssc_prof<ssc(jj).NBssc
        ssc(jj).SSC(ssc(jj).Depth>(bed-2))=ssc(jj).NBssc;
    end
end

load('YR_Jan19_Neap.mat')
% % adcp time and index will be reused in the SedFlux calculation
% % only run this section once, the first time you process
% for jj=1:length(adcp)
%     %get the matlab time for each transect and assign a transect number to
%     %each time point
%     adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
%         adcp(jj).hour,adcp(jj).minute,zeros(1,length(adcp(jj).minute)));
%     adcp(jj).transidx=jj*ones(1,length(adcp(jj).time));
% end
% save('YR_Jan19_Neap','adcp') 

% pull out all ctd profiles from this day
ssc=ssc(1:12);

allidx_adcp=[2,4,6,8,10,11,12,13,15,17,19,22];
for jj=1:length(ssc)
    % get the adcp transect
    ssc(jj).transect=allidx_adcp(jj);
end

%% find idealized transect end points (from adcp processing:
% BR_neap_Residuals)
[right.x,right.y]=deg2utm(16.758914, 96.168642);
[left.x,left.y]=deg2utm(16.766158, 96.170355);
for jj=1:length(ssc)   
    % convert the coordinates of your ADCP records to UTM
    [ssc(jj).x, ssc(jj).y,~] = deg2utm(ssc(jj).lat, ssc(jj).long);
    % project the ADCP data onto the idealized transect
    sscproj = proj([left.x - ssc(jj).x, left.y - ssc(jj).y],...
        [left.x- right.x, left.y - right.y]);
    % define new distances along the idealized transect
    ssc(jj).dist =  round(sqrt(sscproj(1).^2 + sscproj(2).^2));
    
end
% now you have indexed CTD profiles, put those into useful grid matching
% the adcp style:
% 0:1000 wide (1 m bin width)
[~,cols]=size(adcp(1).alongComplete);
% 0:0.01:20m deep (change to exceed max depth of adcp data)
interpdepth=0:0.05:20;

%%
for jj=1:length(ssc)
    ssc(jj).interpalong=NaN([length(interpdepth),cols]);
    % add the ssc profile to the correct location along the transect and
    % interpolate it onto a consistent depth matrix
    ssc(jj).interpalong(:,ssc(jj).dist)=interp1(ssc(jj).Depth,ssc(jj).SSC,interpdepth);
    
    % fill the surface and bed with avg of nearest value
    fillsurface=find(interpdepth<min(ssc(jj).Depth));
    ssc(jj).interpalong(fillsurface,ssc(jj).dist)=nanmean(...
        ssc(jj).interpalong(fillsurface(end):fillsurface(end)+3,ssc(jj).dist));
    fillbed=find(interpdepth>max(ssc(jj).Depth));
    ssc(jj).interpalong(fillbed,ssc(jj).dist)=nanmean(...
        ssc(jj).interpalong(fillbed(1)-3:fillbed(1),ssc(jj).dist));
    ssc(jj).adcpbeddepth=adcp(ssc(jj).transect).interpdepths;
    

clear fill*
end
clear *adcp iC iA T

%% Put the ctd profiles into sigma coordinates in a new matrix:
% 100pt vector between 0 and 1 for the sigma depths
% calc the sigma depths for each column (divide depth by max depth)
% interpolate onto the 100pt matrix
sigmadepth=linspace(0,1,50);

for jj=1:length(ssc)
    % fill to bed with nearest value
    ssc(jj).interpalong=fillmissing(ssc(jj).interpalong,'nearest',1);
    
    % divide depth by max depth:
    ssc(jj).sigmaDepth=interpdepth'./ssc(jj).adcpbeddepth;
    % make a nan matrix that will be the sigma interpolation matrix
    ssc(jj).sigmaAlongComplete=NaN([length(sigmadepth),1001]);
    
    for nn=ssc(jj).dist
        ssc(jj).sigmaAlongComplete(:,nn)=interp1(ssc(jj).sigmaDepth(:,nn)...
            ,ssc(jj).interpalong(:,nn),sigmadepth);
    end
end


%% next interpolate ssc profiles across transects:

for jj=1:length(ssc)
    % very simple QC/data clean up: ssc should be positive
    ssc(jj).sigmaAlongComplete(ssc(jj).sigmaAlongComplete<=0)=NaN;
    
    % fill linearly across river
    ssc(jj).sigmaAlongComplete=fillmissing(ssc(jj).sigmaAlongComplete,'nearest',2);
    
    figure;
    pcolor(0:1000,sigmadepth,ssc(jj).sigmaAlongComplete),shading flat, colorbar
end


save('YR_Jan19_Neap_SSC','ssc')