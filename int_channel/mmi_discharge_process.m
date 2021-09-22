%% make a struc that is just from interior island surveys
clear all,close all,clc

cd C:\GLOVER\data\myanmar\mmi_velocity\BR_Int

% Folder containing only text files of each transect
FileList=dir('*ASC.TXT');

for n = 1:length(FileList)
    disp(['processing ' num2str(n)]);
    
    adcp(n) = adcp_ascii(FileList(n).name);

end

for jj=1:length(adcp)
        adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
        adcp(jj).hour,adcp(jj).minute,adcp(jj).second);
end
cd C:\GLOVER\output\myanmar\mmi_discharge
save('br_08mar2018','adcp')
%%
clear all,close all,clc

cd C:\GLOVER\data\myanmar\mmi_velocity\EastChannel_Mar18
% Folder containing only text files of each transect
FileList=dir('*ASC.TXT');

for n = 1:length(FileList)
    disp(['processing ' num2str(n)]);
    
    adcp(n) = adcp_ascii(FileList(n).name);
end

for jj=1:length(adcp)
    
        adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
        adcp(jj).hour,adcp(jj).minute,adcp(jj).second);
end

cd C:\GLOVER\output\myanmar\mmi_discharge
save('br_07mar2018','adcp')

%%
clear all,close all,clc
cd C:\GLOVER\output\myanmar\mmi_discharge

load('br_08mar2018')
int = adcp(1);
load('br_07mar2018')

int = [int;adcp(8);adcp(9)];
clear adcp
adcp = int;

for jj=1:3
    adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
        adcp(jj).hour,adcp(jj).minute,adcp(jj).second);
    t(jj)=adcp(jj).time(1);
end
[t,idx]=sort(t,'ascend');
adcp=adcp(idx);

% Fill ADCP data for top and bottom of water column
for i=1:length(adcp)
    bottom2add=[];
    cols=[];
    columns=[];
    measuredBins=[];
    n=[];
    rows=[];
    temprow=[];
   
    % find the avg bed depth
    adcp(i).depth=nanmean(adcp(i).depth);
    adcp(i).depth(adcp(i).depth<0)=NaN;
    adcp(i).depth=fillmissing(adcp(i).depth,'nearest'); 
    
    % determine size of interpreted latrix
    [rows,cols]=size(adcp(i).spd);
    
    % define a new z matrix to allow for upper water column
    adcp(i).zComplete=nan(rows+2,cols);
    adcp(i).zComplete(3:rows+2,:)=repmat(adcp(i).z(:,1),1,length(adcp(i).zComplete),1);
    adcp(i).zComplete(1,:)=adcp(i).z(1)-1;
    adcp(i).zComplete(2,:)=adcp(i).z(1)-0.5;
    
    % fill in the upper water column with the mean velocity from the upper
    % two measured bins
    adcp(i).spdComplete=nan(rows+2,cols);
    adcp(i).spdComplete(3:rows+2,:)=adcp(i).spd;
%     adcp(i).spdComplete(1,:)=nanmean(adcp(i).spdComplete(3:4,:));
    adcp(i).spdComplete(2,:)=nanmean(adcp(i).spdComplete(3:4,:));
    
     % fill the bottom bins           
     [~,temprow]=max(~isnan(flipud(adcp(i).spdComplete)), [], 1); %find first bin with nan data as the max of measured values
     temprow(temprow==1)=NaN; %if this is the first bin then ignore that column (empty)
     measuredBins=rows+2-temprow+1; %for the along complete, the row index for this value will be the number of rows+2 for surface minus the 
     bottom2add=(round(adcp(i).depth)-measuredBins./2)./(adcp(i).binsize./100);
     [~,columns]=find(measuredBins>1);

     % New method (HEG, 8/30/18), fill with 1/2 and 1/4 avg for every 2 bins
    for n=columns
        % fill first 2 missing bins with avg of deepest 2 measured
        fillavg=nanmean(adcp(i).spdComplete(measuredBins(n)-2:measuredBins(n),n));
        adcp(i).spdComplete(measuredBins(n)+1:measuredBins(n)+2,n)=fillavg;
        if bottom2add(n)<=4 %last 2 of 4 missing bins is half of the above avg
            adcp(i).spdComplete(measuredBins(n)+3:measuredBins(n)+bottom2add(n),n)=fillavg/2;
        else %if more than 3 bins, first 2 are avg of above 2 and rest are half of that avg
            adcp(i).spdComplete(measuredBins(n)+3:measuredBins(n)+4,n)=fillavg/2;
            adcp(i).spdComplete(measuredBins(n)+5:measuredBins(n)+bottom2add(n),n)=fillavg/4;     
        end
    end
%     % add bottom and top row to smooth then remove those same rows
%     spd = adcp(i).spdComplete;
%     for n=columns
%         spd(measuredBins(n)+bottom2add(n)+1,n) =...
%             spd(measuredBins(n)+bottom2add(n),n);
%     end
%     spd = [spd(1,:);spd];
%     K = (1/25)*ones(5);
%     spd = conv2(spd,K,'same');
%     
%     adcp(i).spdComplete = spd(2:end,:);
    
%     figure;pcolor(adcp(i).time,adcp(i).zComplete,adcp(i).spdComplete),shading flat,axis ij,colorbar

end

save('br_int_march2018','adcp')


%% add in CTD profiles of turbidity and salinity:

%% pull out ctd profiles from Int channel survey March 2018
clear all,close all,clc
cd C:\GLOVER\output\myanmar

load('AyeMar18_CTD_all.mat')

for jj=1:length(CombinedProfiles)
    
    t(jj)=floor(CombinedProfiles(jj).time(1));
end
[t,idx]=sort(t,'ascend');
CombinedProfiles=CombinedProfiles(idx);
ctd=CombinedProfiles(t==datenum('8-March-2018'));

save('C:\GLOVER\output\myanmar\mmi_discharge\ctdprof_8Mar18.mat','ctd')

%% pull out ctd profiles from Int channel survey Sep 2017
clear all,close all,clc
cd C:\GLOVER\output\myanmar

load('AyeSept17_CTD_all.mat')

for jj=1:length(CombinedProfiles)
    
    t(jj)=floor(CombinedProfiles(jj).time(1));
end
[t,idx]=sort(t,'ascend');
CombinedProfiles=CombinedProfiles(idx);
ctd=CombinedProfiles(t==datenum('13-Sept-2017'));

save('C:\GLOVER\output\myanmar\mmi_discharge\ctdprof_13Sep17.mat','ctd')

%% add dist along chnl to ctd profile struct
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('mmi_discharge\ctdprof_13Sep17.mat')
load('mmi_discharge\channel_trace.mat')

%figure out where ctd profs are along channel

for jj=1:length(ctd)
    [x,y,~] = deg2utm(ctd(jj).lat, ctd(jj).long);
    dst = sqrt(abs(ch.x-x).^2 + abs(ch.y-y).^2);
    [~,idx]=min(dst);
    ctd(jj).dist_in=ch.dist_in(idx);
    ctd(jj).dist=ch.dist(idx);
end
save('C:\GLOVER\output\myanmar\mmi_discharge\ctdprof_13Sep17.mat','ctd')
%% add dist along chnl to ctd profile struct
clear all,close all,clc
cd C:\GLOVER\output\myanmar
load('mmi_discharge\ctdprof_8Mar18.mat')
load('mmi_discharge\channel_trace.mat')

%figure out where ctd profs are along channel

for jj=1:length(ctd)
    [x,y,~] = deg2utm(ctd(jj).lat, ctd(jj).long);
    dst = sqrt(abs(ch.x-x).^2 + abs(ch.y-y).^2);
    [~,idx]=min(dst);
    ctd(jj).dist_in=ch.dist_in(idx);
    ctd(jj).dist=ch.dist(idx);
end
save('C:\GLOVER\output\myanmar\mmi_discharge\ctdprof_8Mar18.mat','ctd')