clear all,close all,clc

% compare elevation to water level
load('C:\GLOVER\data\myanmar\SurveyPts2019.mat')
cd C:\GLOVER\output\myanmar
load('aqd\aqd_sep19_hc.mat')
idx = find(survey.ID==1 & survey.type==0);

% pull out all data for croc stn
flds = fieldnames(survey);
for jj=1:7
    croc.(flds{jj})=survey.(flds{jj})(idx);
end
% BR34-36 at '9/28/2019 8:20 AM';'9/28/2019 8:20 AM';'9/28/2019 8:10 AM'};
wl_time=datenum('9/28/2019 8:10 AM');
wl_elev = croc.Zcorr(end);
% wl_dist = croc.dist(end);

right.x = max(croc.Lat);
right.y = max(croc.Lon);
left.x = min(croc.Lat);
left.y = min(croc.Lon);

for jj=1:length(croc.Lat)
    p = proj([left.x-croc.Lat(jj), left.y-croc.Lon(jj)],...
        [left.x-right.x, left.y-right.y]);
    croc.dist(jj) =  sqrt(p(1).^2 + p(2).^2);
end
[croc.dist,srt]=sort(croc.dist,'ascend');

for jj=1:7
    croc.(flds{jj})=croc.(flds{jj})(srt);
end
croc.dist = round(croc.dist,1);
dist=0:0.1:120.2;
[~,ia,ib] = intersect(croc.dist,dist);
elev = NaN(size(dist));
elev(ib) = croc.Zcorr(ia);
elev(elev<-43.5)=NaN;
elev(dist==100)=elev(end)-aqd.depth(aqd.time==wl_time);
elev(dist==90)=elev(end)-5;
elev(dist==70)=elev(end)-3.7;
elev(dist==60)=elev(end)-3.7;

elev = fillmissing(elev,'linear');

ht = elev(end)+(max(aqd.depth) - aqd.depth(aqd.time==wl_time));
lt = elev(end)+(min(aqd.depth) - aqd.depth(aqd.time==wl_time));

figure;
plot(croc.dist,croc.Zcorr,'o'),hold on
plot(dist,elev,'-')
refline(0,ht),refline(0,lt)


figure;plot(aqd.depth)
%%
clear all,close all,clc

F=dir('BR_*_FluxDecomp3.mat');

for jj=1:length(F)
    load(F(jj).name)
    figure;plot(fluxdecomp.A)
    F(jj).dA = max(fluxdecomp.A)-min(fluxdecomp.A);
end