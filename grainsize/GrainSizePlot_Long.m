%% Make Grain Size figure for manuscript:
% plot the low-flow
clear all,close all,clc
rvr={'YR';'BR';'PR'};
for kk=1%:3
load([rvr{kk},'_GrainSize_Sept17.mat'])

% project the ctd locations to the thalweg to get distance along estuary
lat=vertcat(D.lat);
lon=vertcat(D.long);
dist = round(AYEdistance(lat,lon,kk)/1000);
% make new phi scale to interpolate onto
phi = -log(0.37/1000):-0.01:-log2(2);
for jj=1:length(D)
    
    % convert to phi and linearly interp to improve resolution:
    D(jj).data_phi(:,1)=phi;  % put in new phi bins  
    D(jj).data(:,1)=-log2(D(jj).data(:,1)/1000);% convert orig bins to phi
    D(jj).data_phi(:,2)=interp1(D(jj).data(:,1),D(jj).data(:,2),phi)'; %interp
    D(jj).data_phi(:,2) = ceil(D(jj).data_phi(:,2)./...
        sum(D(jj).data_phi(:,2),'omitnan').*10000);% convert to counts out of 10000
    % remove zeros
    D(jj).data_phi(D(jj).data_phi(:,2)==0 | isnan(D(jj).data_phi(:,2)),:)=[];
    % replicate phi by num of counts
    D(jj).counts = repelem(D(jj).data_phi(:,1),D(jj).data_phi(:,2));

    mm(jj)=length(D(jj).counts); % find length of new counts vector
    
    % D50:
    d50(jj) = D(jj).phi.median;
end

% sort data and distance by distance
[dist,idx] = sort(dist,'ascend');
d50=d50(idx); D=D(idx); mm=mm(idx);

H = NaN(max(mm),length(dist));
for jj=1:length(D)
    H(1:mm(jj),jj)=D(jj).counts;
end

dist_unq=unique(dist);
%dist_unq=unique(round(dist,-3));
for jj=1:length(dist_unq)
    Hunq(:,jj)=nanmean(H(:,dist==dist_unq(jj)),2);
    d50unq(jj)=nanmean(d50(dist==dist_unq(jj)));
end

figure(1)
subplot(3,1,kk)
[~,~,~,~,~]=violin(Hunq,'x',dist_unq-0.1,'MED',d50unq,...
    'facecolor',[0.2 0.2 0.7],'facealpha',0.9);
hold on
% figure(2)
% subplot(3,1,kk)
% [~,~,~,~,~]=violin(H,'x',dist,'MED',d50,'facecolor',[0.8 0.8 0.8]);
clearvars -except kk rvr
end

%% plot high flow
clear all
rvr={'YR';'BR';'PR'};
for kk=1:3
load([rvr{kk},'_GrainSize_Mar18.mat'])

% project the ctd locations to the thalweg to get distance along estuary
lat=vertcat(D.lat);
lon=vertcat(D.long);
dist = round(AYEdistance(lat,lon,kk)/1000);

% make new phi scale to interpolate onto
phi = -log(0.37/1000):-0.01:-log2(2);
for jj=1:length(D)
  
    % convert to phi and linearly interp to improve resolution:
    D(jj).data_phi(:,1)=phi;  % put in new phi bins  
    D(jj).data(:,1)=-log2(D(jj).data(:,1)/1000);% convert orig bins to phi
    D(jj).data_phi(:,2)=interp1(D(jj).data(:,1),D(jj).data(:,2),phi)'; %interp
    D(jj).data_phi(:,2) = ceil(D(jj).data_phi(:,2)./...
        sum(D(jj).data_phi(:,2),'omitnan').*10000);% convert to counts out of 10000
    % remove zeros
    D(jj).data_phi(D(jj).data_phi(:,2)==0 | isnan(D(jj).data_phi(:,2)),:)=[];
    % replicate phi by num of counts
    D(jj).counts = repelem(D(jj).data_phi(:,1),D(jj).data_phi(:,2));

    mm(jj)=length(D(jj).counts); % find length of new counts vector
    
    % D50:
    d50(jj) = D(jj).phi.median;
end

% sort data and distance by distance
[dist,idx] = sort(dist,'ascend');
d50=d50(idx); D=D(idx); mm=mm(idx);

H = NaN(max(mm),length(dist));
for jj=1:length(D)
    H(1:mm(jj),jj)=D(jj).counts;
end

dist_unq=unique(dist);
for jj=1:length(dist_unq)
    Hunq(:,jj)=nanmean(H(:,dist==dist_unq(jj)),2);
    d50unq(jj)=nanmean(d50(dist==dist_unq(jj)));
end

figure(1)
subplot(3,1,kk)
[~,~,~,~,~]=violin(Hunq,'x',dist_unq+0.1,'MED',d50unq,...
    'facecolor',[0.3 0.3 0.3]);
ylim([-1 9])
clearvars -except kk rvr

end

subplot(3,1,1),xlim([5 45])
subplot(3,1,2),xlim([5 45])
subplot(3,1,3),xlim([57 100])