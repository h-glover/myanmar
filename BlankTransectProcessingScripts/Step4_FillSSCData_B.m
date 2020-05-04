% Use this code to convert individual CTD casts from the compiled mat into
% a series of interpolated transects associated with an adcp transect. Use
% B when you have 1 cast per transect. This code first compiles the casts
% by transect, then converts to sigma coordinates, then interpolates.
% 1. Associate cast with its transect
% 2. put casts onto an individual transect matrix: X rows deep based on deepest cast and y columns wide based on standard river width.
% 3. Project casts onto transect based on fixed middle pt along transect
% 4. Pad this matrix with values to seabed
% 5. Convert to sigma coordinates
% 6. Interpolate across transects

clear all,close all,clc
load('AyeJan19_CTD_all.mat')

ssc=profiles; clear profiles
for jj=1:length(ssc)
    ssc(jj).profileidx=jj; %label 1:length(ssc) for later use
    ssc(jj).timeidx=ssc(jj).time(1);%time cast starts, for matching adcp transect
end
% concatenate the profile indexs and times
allidx_ssc=horzcat(ssc.profileidx);
alltime_ssc=vertcat(ssc.timeidx);
alltime_ssc=datevec(alltime_ssc);alltime_ssc(:,6)=0;
alltime_ssc=datenum(alltime_ssc);


load('adcp.mat')
% % adcp time and index will be reused in the SedFlux calculation
% % only need to run this section once, the first time you process
% for jj=1:length(adcp)
%     %get the matlab time for each transect and assign a transect number to
%     %each time point
%     adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
%         adcp(jj).hour,adcp(jj).minute,zeros(1,length(adcp(jj).minute)));
%     adcp(jj).transidx=jj*ones(1,length(adcp(jj).time));
% end
% save('YR_Jan19_Neap','adcp') 

% concatenate transect index into 1 vector 
allidx_adcp=horzcat(adcp.transidx);
alltime_adcp=horzcat(adcp.time);
% compare adcp time and ctd time (minutes removed)
[~,iA,iC]=intersect(alltime_adcp,alltime_ssc);
% pull out all ctd profiles from this day
ssc_trans=allidx_ssc(iC);
allidx_adcp=allidx_adcp(iA);
ssc=ssc(iC);
for jj=1:length(ssc)
    % get the adcp transect
    ssc(jj).transect=allidx_adcp(jj);
end

% find idealized transect end points (from adcp processing:
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

for jj=1:length(ssc)
    ssc(jj).interpalong=NaN([length(interpdepth),cols]);
    % intersect ctd cast depth with the interpolated depth for the big matrix
    [~,i1,i2]=intersect(round(ssc(jj).Depth,2),interpdepth);
    % fill in matrix column using distance across transect and interpolated
    % depth index
    ssc(jj).interpalong(i2,ssc(jj).dist)=ssc(jj).SSC(i1);
    % find the min and max depth of the profile and interpolate within that range
    fillidx=find(interpdepth>=min(ssc(jj).Depth) & interpdepth<=max(ssc(jj).Depth));
    ssc(jj).interpalong(fillidx,ssc(jj).dist)=fillmissing(...
        ssc(jj).interpalong(fillidx,ssc(jj).dist),'linear');
    ssc(jj).interpalong(1:fillidx(1),ssc(jj).dist)=nanmean(...
        ssc(jj).interpalong(fillidx(1):fillidx(1)+5,ssc(jj).dist));
    ssc(jj).adcpbeddepth=round(adcp(ssc(jj).transect).interpdepths,2);
    clear i1 i2
end
clear *adcp iC iA T

%% Put the ctd profiles into sigma coordinates in a new matrix:
% 100pt vector between 0 and 1 for the sigma depths
% calc the sigma depths for each column (divide depth by max depth)
% interpolate onto the 100pt matrix
sigmadepth=linspace(0,1,100);

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
    pcolor(0:1000,sigmadepth,ssc(jj).sigmaAlongComplete)
    shading flat,axis ij,colorbar
end



save('','ssc')