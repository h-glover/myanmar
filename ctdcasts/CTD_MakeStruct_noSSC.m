
%%
% Files you'll need:
% - RBR CTD data exported from Ruskin via the "Legacy" export (Export Dataset>Legacy>Matlab MAT)
% - MAT file containing the fllowing structures:
%       - GPSStruct (structure containing all GPS points from cruise) with these fields:
%           - Datetime, Lat, Long, WP, and Time (Datetime in Matlab
%            dateimt format, and Time in serial date format).
%       SSCStruct (structure containing info from filtered SSC samples)
%       with these fields:
%           - WP, SSC, Lat, Long, Datetime, and Time Datetime in Matlab
%            dateimt format, and Time in serial date format). 

% Make sure the .mat file exported by Ruskin is in your path, and enter the
% file name below:
clear all,close all,clc
load('018042_20201011_1611.mat')

for jj=1:2
    RBR.channelnames{jj+3}=sprintf('Voltage %d',jj);
end

fid = fopen('GPS_CTD.csv');
T=textscan(fid,'%s %s %s %s %f %f','Delimiter',',');
fclose(fid);

GPSStruct.Lat=T{5}(9:11);
GPSStruct.Long=T{6}(9:11);
GPSStruct.Wp=T{4}(9:11);
GPSStruct.Name=T{3}(9:11);


%%


% Make sure the .mat file containing GPSStruct and SSCStruct is in your path, and enter the
% file name below:

% Plot depth data 
fig = figure(1);
set(gcf,'color','w')
clf
ax=tight_subplot(1,1,.1,.1,.05);
axes(ax(1))
plot(RBR.data(:,3),'k','linewidth',2),axis ij,hold on
xlim([1,length(RBR.data)])
hold on


% numplots=4;
% ax=tight_subplot(numplots,1,.1,.1,.1);
% c=0;
% for n=1:numplots
% axes(ax(n))
% plot(RBR.data(1+(c*1/numplots)*length(RBR.data):((c+1)*1/numplots)*length(RBR.data),3),'k','linewidth',2),axis ij,hold on
% c=c+1;
% end
% plot(RBR.data(:,9),'*k')
%% Selecting points

numcasts=input('Number of Casts?');% Count the number of downcasts, and enter it when prompted
%21 profiles
% pick the rough location of each cast
disp('Click the approximate location of each cast')

% collect rough cast location info from user
for kk=1:numcasts
    disp(kk)
    [centers(1,kk),~] = ginput(1); disp(centers(1,kk))
end



% zoom plot, and allow user to click the beginning and end of each cast
for jj=1:numcasts
    disp(jj)
    
    figure(1)
         xlim([centers(jj)-400,centers(jj)+200])
 
    s=scatter(centers(jj),RBR.data(round(centers(jj)),3),100,'r');
%     pause(1)
    [pickpt(1,jj),~] = ginput(1); disp(pickpt(1,jj))%pick start surf
%     pause(5)
    [pickpt(2,jj),~] = ginput(1); disp(pickpt(2,jj))%pick end bottom
    s=scatter(centers(jj),RBR.data(round(centers(jj)),3),100,'w');
 

end

% store the user's selections for beginning and end of casts.
pickpt=round(pickpt);
save('pickpt','pickpt')


%% Build the data structure
load('pickpt2.mat') %load the points if you previously picked them
varnames=RBR.channelnames';
varnames=strrep(varnames,' ','');
RBR.datenum=datenum(RBR.sampletimes,'dd/mm/yyyy HH:MM:SS.FFF');
RBR.datenum=datevec(RBR.datenum); RBR.datenum(:,4)=RBR.datenum(:,4)+6;
RBR.datenum(:,5)=RBR.datenum(:,5)+30; RBR.datenum=datenum(RBR.datenum);
RBR.datetime=datetime(RBR.datenum,'ConvertFrom','datenum');
DepthColumn=find(strcmp(varnames, 'Depth'));
RBR.data(:,DepthColumn)=RBR.data(:,DepthColumn)-min(RBR.data(:,DepthColumn));

[~,numcasts]=size(pickpt);

% Build the structure
for jj=1:length(pickpt)
    for ii=1:length(varnames)
        profiles(jj).(varnames{ii})=RBR.data(pickpt(1,jj):pickpt(2,jj),ii); % add all the data from the RBR mat file
    end
    profiles(jj).time=RBR.datenum(pickpt(1,jj):pickpt(2,jj)); % add time info
    profiles(jj).datetime=RBR.datetime(pickpt(1,jj):pickpt(2,jj)); % add a datetime field
    profiles(jj).lat=GPSStruct.Lat(jj);
    profiles(jj).long=GPSStruct.Long(jj);
    profiles(jj).wp=GPSStruct.Wp(jj);
    profiles(jj).name=GPSStruct.Name(jj);
end

[~,filename,~] = fileparts(RBR.datasetfilename);
filename=strcat(filename,'_profiles.mat');
save(filename,'profiles','pickpt');
%% Plot a figure showing the turbidity profiles

for jj=1:length(profiles)
    figure;
    subplot(121)
    plot(profiles(jj).SSC,profiles(jj).Depth)
    axis ij
    subplot(122)
    plot(profiles(jj).Salinity,profiles(jj).Depth)
    axis ij
end