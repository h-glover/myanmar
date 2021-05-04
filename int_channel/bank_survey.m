% clear all,close all,clc
% 
% cd C:\GLOVER\output\myanmar\mmi_discharge
% load('channel_trace.mat')
% fid = fopen('banksurvey.csv');
% t = textscan(fid,'%f %f %f %s %f','Delimiter',',');
% 
% surv.lat=t{1};
% surv.long=t{2};
% surv.wpt=t{3};
% surv.side=zeros(56,1);
% surv.side(strcmp(t{4},'w'))=1;
% surv.width=t{5};
% surv.key={'w=1','e=0'};
% surv.dist=NaN(56,1);
% for jj=1:56
%     [x,y,~] = deg2utm(surv.lat(jj), surv.long(jj));
%     dst = sqrt(abs(ch.x-x).^2 + abs(ch.y-y).^2);
%     [~,idx]=min(dst);
%     surv.dist_in(jj)=ch.dist_in(idx);
%     surv.dist(jj)=ch.dist(idx);
% end
% 
% save('BankSurvey','surv')

%%
% small: 0-5 m
% medium: 5-10 m
% large: 10-15 m
clear all,close all,clc
cd C:\GLOVER\output\myanmar\mmi_discharge\
load('channel_trace.mat')
load('BankSurvey.mat')
load('br_int_march2018.mat'),adcp=adcp(3);

adcp.depth = movmean(adcp.depth,7);

surv.width(:,2:4)=0; %col1 = small, col2 = med, col3 = large
surv.width(surv.width(:,1)<=5,2) = 1;
surv.width(surv.width(:,1)>5 & surv.width(:,1)<=10,3) = 1;
surv.width(surv.width(:,1)>10,4) = 1;

% calc channels per km
km = 0:1000:14000;
for jj=2:length(km)
    % km_ch(jj-1) = length(find(surv.dist>km(jj-1) & surv.dist<km(jj)));
    ch_dist(jj-1,1:3) = sum(surv.width(surv.dist>km(jj-1) & surv.dist<km(jj),2:4));
end
% km = km(2:end)/1000;
km = (500:1000:km(end))/1000;

figure;
subplot(211)
yyaxis left,scatter(surv.dist/1000,surv.width(:,1),[],surv.width(:,1),'filled')
caxis([0 10]),colorbar,colormap(cmocean('turbid'))
ylabel('channel size')
yyaxis right,plot(adcp.elapdist/1000,adcp.depth,'k')
ax=gca;ax.YColor = 'k';axis ij,ylabel('channel depth')
xlim([0 14]),xlabel('km from south entrance')

subplot(212)
yyaxis left
for jj=1:3
    plot(km,ch_dist(:,jj)),hold on
end
ylabel('channels/km')
yyaxis right,plot(adcp.elapdist/1000,adcp.depth,'k')
ax=gca;ax.YColor = 'k';axis ij,ylabel('channel depth')
xlim([0 14]),xlabel('km from south entrance')
legend({'small','med','large','depth'})

figure;
plot(ch.lon,ch.lat,'k-'),hold on
scatter(surv.long,surv.lat,surv.width(:,1)*10,surv.width(:,1),'filled')
colorbar,colormap(cmocean('turbid')),caxis([0 10])