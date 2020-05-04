load YR_Jan19_Neap.mat

% fix the depths for the missing data sections
C=[1/9];
C=repmat(C,[1,9]);
for i=1:length(adcp)
    % fix the adcp depths by smoothing and filling in the gaps with the avg depth
    adcp(i).interpdepths=conv(adcp(i).interpdepths,C,'same');
    adcp(i).interpdepths(300)=18;
    adcp(i).interpdepths(30:800)=fillmissing(...
        adcp(i).interpdepths(30:800),'linear');
end

for i=1:length(adcp)
    bottom2add=[];
    cols=[];
    columns=[];
    measuredBins=[];
    n=[];
    rows=[];
    temprow=[];
    
    % determine size of interpreted matrix
    [rows,cols]=size(adcp(i).interpalong);
    
    % define a new z matrix to allow for upper water column
    adcp(i).zComplete=nan(rows+2,cols);
    adcp(i).zComplete(3:rows+2,:)=repmat(adcp(i).z(:,1),1,length(adcp(i).zComplete),1);
    adcp(i).zComplete(1,:)=adcp(i).z(1)-1;
    adcp(i).zComplete(2,:)=adcp(i).z(1)-0.5;
    
    % fill in the upper water column with the mean velocity from the upper
    % two measured bins
    adcp(i).alongComplete=nan(rows+2,cols);
    adcp(i).alongComplete(3:rows+2,:)=adcp(i).interpalong;
    adcp(i).alongComplete(1,:)=mean(adcp(i).interpalong(3:4,:));
    adcp(i).alongComplete(2,:)=mean(adcp(i).interpalong(3:4,:));
    
    % fill the large gaps in this Yangon record by linear interpolation
    D=round(adcp(i).uniquedist(end));
    absmax=max(max(abs(adcp(i).alongComplete)));
    L=[];
    L=isnan(nanmean(adcp(i).alongComplete(1:30,1:D)));
    L(1:20)=0;
    adcp(i).fillgap=fillmissing(adcp(i).alongComplete,'linear',2);
    adcp(i).alongComplete(1:30,L)=adcp(i).fillgap(1:30,L);
    adcp(i).alongComplete(abs(adcp(i).alongComplete)>absmax)=NaN;
    
    % fill the bottom bins
    [~,temprow]=max(~isnan(flipud(adcp(i).alongComplete)), [], 1); %find first bin with nan data as the max of measured values
    temprow(temprow==1)=NaN; %if this is the first bin then ignore that column (empty)
    measuredBins=rows+2-temprow+1; %for the along complete, the row index for this value will be the number of rows+2 for surface minus the
    bottom2add=(round(adcp(i).interpdepths)-measuredBins./2)./(adcp(i).binsize./100);
    [~,columns]=find(measuredBins>1);
    
    % %      % previous method: all bins filled with constant value
    % %     for n=columns
    % %         adcp(i).alongComplete(measuredBins(n):measuredBins(n)+bottom2add(n),n)=...
    % %             nanmean(adcp(i).alongComplete(measuredBins(n)-2:measuredBins(n),n));
    % %     end
    % New method (HEG, 8/30/18), fill with 1/2 and 1/4 avg for every 2 bins
    for n=columns
        % fill first 2 missing bins with avg of deepest 2 measured
        fillavg=nanmean(adcp(i).alongComplete(measuredBins(n)-2:measuredBins(n),n));
        adcp(i).alongComplete(measuredBins(n)+1:measuredBins(n)+2,n)=fillavg;
        if bottom2add(n)<=4 %last 2 of 4 missing bins is half of the above avg
            adcp(i).alongComplete(measuredBins(n)+3:measuredBins(n)+bottom2add(n),n)=fillavg/2;
        elseif bottom2add(n)>4 %if more than 3 bins, first 2 are avg of above 2 and rest are half of that avg
            adcp(i).alongComplete(measuredBins(n)+3:measuredBins(n)+4,n)=fillavg/2;
            adcp(i).alongComplete(measuredBins(n)+5:measuredBins(n)+bottom2add(n),n)=fillavg/4;
        end
    end
    
    % check the interpolation...
    figure;
    pcolor(0:1000,adcp(i).zComplete(:,1),adcp(i).alongComplete)
    shading flat,colorbar
end

save YR_Jan19_Neap.mat adcp