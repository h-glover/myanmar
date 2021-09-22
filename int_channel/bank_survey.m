%% load external channel data
clear all,close all,clc

cd C:\GLOVER\data\myanmar
fid = fopen('BankSurveyExt_28Sept19.csv');
t = textscan(fid,'%f %f %f','Delimiter',',','HeaderLines',1);
fclose(fid);

% load('C:\GLOVER\output\myanmar\mmi_discharge\channel_trace.mat')

surv.lat=t{2};
surv.long=t{3};
surv.width=t{1};
surv.dist=NaN(38,1);

st.x = min(surv.lat);
st.y = min(surv.long);
[st.x,st.y,~] = deg2utm(st.x,st.y);

for jj=1:38
    [x,y,~] = deg2utm(surv.lat(jj), surv.long(jj));
    surv.dist(jj) = sqrt(abs(st.x-x).^2 + abs(st.y-y).^2);
end
[surv.dist,idx]=sort(surv.dist,'ascend');
surv.lat=surv.lat(idx);
surv.long=surv.long(idx);
surv.width=surv.width(idx);
surv.width(surv.width==50)=0;

[num_width,ch_width] = groupcounts(surv.width);

tot_ch_width = sum(ch_width.*(61*num_width/(max(surv.dist)/1000)));
totch = (sum(num_width)/(max(surv.dist)/1000))*61;

% % divide into channel size categorie
% surv.width(:,2:4)=0; %col1 = small, col2 = med, col3 = large
% surv.width(surv.width(:,1)<=5,2) = 1;
% surv.width(surv.width(:,1)>5 & surv.width(:,1)<=10,3) = 1;
% surv.width(surv.width(:,1)>5,4) = 1;
% 
% avg_chls = [length(surv.lat),sum(surv.width(:,2:4))]/(max(surv.dist)/1000);



%% load internal channel data
% clear all,close all,clc
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

% %% plot grain size in interior channel during HF and LF
% clear all,close all,clc
% cd C:\GLOVER\output\myanmar
% load('mmi_discharge\channel_trace.mat')
% load('grainsize\BR_GrainSize_Sept17.mat'), hf = D([16:17,49,50,54:57]);
% load('grainsize\BR_GrainSize_Mar18.mat'), lf = D([47:58]);
% 
% lat = vertcat(lf.lat);
% lon = vertcat(lf.long);
% 
% for jj=1:length(hf)
%     [x,y,~] = deg2utm(hf(jj).lat, hf(jj).long);
%     dst = sqrt(abs(ch.x-x).^2 + abs(ch.y-y).^2);
%     [~,idx]=min(dst);
%     hf(jj).dist_in=ch.dist_in(idx);
%     hf(jj).dist=ch.dist(idx);
% end
% for jj=1:length(lf)
%     [x,y,~] = deg2utm(lf(jj).lat, lf(jj).long);
%     dst = sqrt(abs(ch.x-x).^2 + abs(ch.y-y).^2);
%     [~,idx]=min(dst);
%     lf(jj).dist_in=ch.dist_in(idx);
%     lf(jj).dist=ch.dist(idx);
% end
% 
% d = vertcat(hf.dist);
% [~,idx] = sort(d);
% hf = hf(idx);
% d = vertcat(lf.dist);
% [~,idx] = sort(d);
% lf = lf(idx);
% 
% 
% save('mmi_discharge\BR_GrainSize_Channel.mat','hf','lf')
%% HC CHANNEL MORPHOLOGY FIGURE FOR PAPER
% small: 0-5 m
% medium: 5-10 m
% large: 10-15 m
clear all,close all,clc
cd C:\GLOVER\output\myanmar\mmi_discharge
load('BR_GrainSize_Channel.mat')
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

% figure;
% plot(ch.lon,ch.lat,'k-'),hold on
% scatter(surv.long,surv.lat,surv.width(:,1)*10,surv.width(:,1),'filled')
% colorbar,colormap(cmocean('turbid')),caxis([0 10])

% smooth depth:
adcp.depth = movmean(adcp.depth,9);

figure;
subplot(311)
yyaxis left,scatter(surv.dist/1000,surv.width(:,1),[],surv.width(:,1),'filled')
caxis([0 10]),colorbar,colormap(cmocean('turbid'))
ylabel('channel size')
yyaxis right,plot(adcp.elapdist/1000,adcp.depth,'k')
ax=gca;ax.YColor = 'k';axis ij,ylabel('channel depth')
xlim([0 14]),xlabel('km from south entrance')

subplot(312)
yyaxis left
for jj=1:3
    plot(km,ch_dist(:,jj)),hold on
end
ylabel('channels/km')
yyaxis right,plot(adcp.elapdist/1000,adcp.depth,'k')
ax=gca;ax.YColor = 'k';axis ij,ylabel('channel depth')
xlim([0 14]),xlabel('km from south entrance')
legend({'<5 m','5-10 m','>10','depth'})
% small: 0-5 m
% medium: 5-10 m
% large: 10-15 m

dist_hf = vertcat(hf.dist)/1000;
dist_lf = vertcat(lf.dist)/1000;
med_hf = vertcat(hf.median);
med_lf = vertcat(lf.median);
frac_hf = vertcat(hf.frac);
frac_lf = vertcat(lf.frac);

subplot(313)
yyaxis left
semilogy(dist_hf,med_hf,'k-'),hold on
semilogy(dist_lf,med_lf,'r-')
ylim([1 200]),ylabel('D_{50} in um')
yyaxis right
plot(dist_hf,frac_hf(:,3),'k:'),hold on
plot(dist_lf,frac_lf(:,3),'r:')
ylim([0 100]),ylabel('% content')
xlim([0 14]),xlabel('km from south entrance')
legend({'high flow','low flow'})











