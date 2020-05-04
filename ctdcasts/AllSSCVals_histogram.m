% Make a histogram of ssc by season
clear all,close all,clc

% vertically compile all CTD turbidity-based SSC and salinity 
% sort rivers based on longitude
load('AyeJan19_CTD_all.mat')
lon=vertcat(profiles.long);
yr_ssc_lf=vertcat(profiles(lon>96).SSC);
pr_ssc_lf=vertcat(profiles(lon<95).SSC);
yr_sal_lf=vertcat(profiles(lon>96).Salinity);
pr_sal_lf=vertcat(profiles(lon<95).Salinity);

load('AyeMar18_CTD_all.mat')
lon=vertcat(CombinedProfiles.long);
yr_ssc_lf=[yr_ssc_lf;vertcat(CombinedProfiles(lon>96).SSCCal)];
br_ssc_lf=vertcat(CombinedProfiles(lon<96 & lon>95).SSCCal);
pr_ssc_lf=[pr_ssc_lf;vertcat(CombinedProfiles(lon<95).SSCCal)];
yr_sal_lf=[yr_sal_lf;vertcat(CombinedProfiles(lon>96).Salinity)];
br_sal_lf=vertcat(CombinedProfiles(lon<96 & lon>95).Salinity);
pr_sal_lf=[pr_sal_lf;vertcat(CombinedProfiles(lon<95).Salinity)];

load('AyeSept17_CTD_all.mat')
lon=vertcat(CombinedProfiles.long);
yr_ssc_hf=vertcat(CombinedProfiles(lon>96).SSCCal);
br_ssc_hf=vertcat(CombinedProfiles(lon<96 & lon>95).SSCCal);
pr_ssc_hf=vertcat(CombinedProfiles(lon<95).SSCCal);
yr_sal_hf=vertcat(CombinedProfiles(lon>96).Salinity);
br_sal_hf=vertcat(CombinedProfiles(lon<96 & lon>95).Salinity);
pr_sal_hf=vertcat(CombinedProfiles(lon<95).Salinity);


% load the water sample SSC values as well
fid=fopen('Aye_SSC_all.csv');
T=textscan(fid,'%s %*s %*f %f %s %f %f %s','HeaderLines',1,'Delimiter',',');
%1=river, 2=ssc, 3=depth, 4=lat, 5=long, 6=season
fclose(fid);

yr_ssc_lf=[yr_ssc_lf;T{2}(T{5}>96 & strcmp('Low',T{6}))];
L(1)=length(yr_ssc_lf);
yr_ssc_hf=[yr_ssc_hf;T{2}(T{5}>96 & strcmp('High',T{6}))];
L(2)=length(yr_ssc_hf);

br_ssc_lf=[br_ssc_lf;T{2}(T{5}<96 & T{5}>95 & strcmp('Low',T{6}))];
L(3)=length(br_ssc_lf);
br_ssc_hf=[br_ssc_hf;T{2}(T{5}<96 & T{5}>95 & strcmp('High',T{6}))];
L(4)=length(br_ssc_hf);

pr_ssc_lf=[pr_ssc_lf;T{2}(T{5}<95 & strcmp('Low',T{6}))];
L(5)=length(pr_ssc_lf);
pr_ssc_hf=[pr_ssc_hf;T{2}(T{5}<95 & strcmp('High',T{6}))];
L(6)=length(pr_ssc_hf);
%%
figure;
xbins=0:10:5000;
histogram(yr_ssc_lf,'Normalization','probability','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','b'),hold on
histogram(yr_ssc_hf,'Normalization','probability','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','b','LineStyle','--')
histogram(br_ssc_lf,'Normalization','probability','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','k')
histogram(br_ssc_hf,'Normalization','probability','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','k','LineStyle','--')
histogram(pr_ssc_lf,'Normalization','probability','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','r')
histogram(pr_ssc_hf,'Normalization','probability','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','r','LineStyle','--')
%%
figure;
xbins=0:4500;
histogram(yr_ssc_lf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','b'),hold on
histogram(yr_ssc_hf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','b','LineStyle','--')
histogram(br_ssc_lf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','k')
histogram(br_ssc_hf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','k','LineStyle','--')
histogram(pr_ssc_lf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','r')
histogram(pr_ssc_hf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','r','LineStyle','--')
xlim([0 4500])
xlabel('SSC (mg/L)')
ylabel('Cumulative Probability')
legend({'Yangon, Low Flow','Yangon, High Flow',...
    'Bogale, Low Flow','Bogale, High Flow',...
    'Pathein, Low Flow','Pathein, High Flow'})
%%
figure;
xbins=0:1000;
histogram(yr_ssc_hf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','b','LineStyle','--'),hold on
histogram(br_ssc_hf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','k','LineStyle','--')
histogram(pr_ssc_hf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','r','LineStyle','--')
xlim([0 1000])
xlabel('SSC (mg/L)')
ylabel('Cumulative Probability')
legend({'Yangon, High Flow','Bogale, High Flow','Pathein, High Flow'})

%
figure;
xbins=0:4500;
histogram(yr_ssc_lf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','b','LineStyle','-'),hold on
histogram(br_ssc_lf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','k','LineStyle','-')
histogram(pr_ssc_lf,'Normalization','cdf','BinEdges',xbins,...
    'DisplayStyle','stairs','EdgeColor','r','LineStyle','-')
xlim([0 4500])
xlabel('SSC (mg/L)')
ylabel('Cumulative Probability')
legend({'Yangon, Low Flow','Bogale, Low Flow','Pathein, Low Flow'})