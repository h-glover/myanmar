% BR neap sediment flux:

% get an ssc profile for each transect: by time, nearest previous
% convert ssc and adcp to sigma coordinates:
% multiply by max depth, interpolate to 1000 point 0:1 range
% transects length stays the same
% multiply ssc by discharge for each transect
clear all,clc

load('BR_Sept17_Neap.mat')
load('BR_Sept17_Neap_SSC.mat')

% get the adcp transect indicies from the ssc transects
sscidx=vertcat(ssc.transect);
% Make a 100x1000 cross section for the sigma format
sigmadepth = linspace(0,1,100).';
L=length(sigmadepth);

for jj=1:length(adcp)
    
    % if there is an ssc matching the current transect, put it in,
    % otherwise, put in the next matching ssc transect
    idx=find(sscidx==jj);
     if isempty(idx)==1
        diff=abs(jj-sscidx);
        [~,idx]=min(diff);
    end
    adcp(jj).ssc=ssc(idx).sigmaAlongComplete;
   
    % put in the CTD locations for future reference
    adcp(jj).CTDlocations=ssc(idx).dist;
    
    % divide adcp depths by the max depth
    adcp(jj).sigmaDepth=adcp(jj).zComplete./adcp(jj).interpdepths;

    % for each column with depth>0, interpolate to the sigma depth
    % vector (100 linearly spaced pts between 0:1)
    for nn=1:length(adcp(jj).interpdepths)
        if isnan(adcp(jj).interpdepths(nn))==1
            adcp(jj).sigmaAlongComplete(:,nn)=NaN([L,1]);
        else
            adcp(jj).sigmaAlongComplete(:,nn)=interp1(...
                adcp(jj).sigmaDepth(:,nn),adcp(jj).alongComplete(:,nn),sigmadepth);
        end       
    end
    % fill the sigma velocity to the surface
    adcp(jj).sigmaAlongComplete=fillmissing(adcp(jj).sigmaAlongComplete,'nearest',1);
    
    % convert velocity to discharge in m3/s: cm to m (/100), m2 to m3 (*binsize)
    binsize=adcp(jj).interpdepths/(L+3); %sigma bin size
    adcp(jj).sigmaAlongComplete=adcp(jj).sigmaAlongComplete./100.*binsize;
    adcp(jj).alongComplete=adcp(jj).alongComplete./100.*0.5; %bin size is 0.5m
    
    % multiply ssc x discharge and convert m3/s to L/s = mg/s
    adcp(jj).sigmaSedFlux=adcp(jj).ssc.*(1000.*adcp(jj).sigmaAlongComplete);
    
    % calculate total discharge for the transects:
    adcp(jj).time=nanmean(adcp(jj).time);
    adcp(jj).MeasuredDis=nansum(nansum(adcp(jj).alongComplete));
    adcp(jj).sigmaMeasuredDis=nansum(nansum(adcp(jj).sigmaAlongComplete));
    adcp(jj).sigmaMeasuredSedFlux=nansum(nansum(adcp(jj).sigmaSedFlux));
    
    % convert sed flux to tonnes/s
    adcp(jj).sigmaMeasuredSedFlux=adcp(jj).sigmaMeasuredSedFlux./(1e9);
end
save('BR_Sept17_Neap_SedFlux','adcp')

%%

