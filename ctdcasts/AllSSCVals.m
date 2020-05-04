% plot all of the CTD turbidity values and all of the SSC sample values on
% a figure like the one Aaron previously made
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


figure;
scatter(yr_sal_hf,yr_ssc_hf,'rd'),hold on
scatter(yr_sal_lf,yr_ssc_lf,'rd')

scatter(pr_sal_hf,pr_ssc_hf,'gx')
scatter(pr_sal_lf,pr_ssc_lf,'gx')

scatter(br_sal_hf,br_ssc_hf,'k*')
scatter(br_sal_lf,br_ssc_lf,'k*')

ax=gca; ax.YScale='Log';ylim([10 10000])

%% load the water sample SSC values as well
fid=fopen('Aye_SSC_all.csv');
T=textscan(fid,'%s %*s %*f %f %s %f %f %s','HeaderLines',1,'Delimiter',',');
%1=river, 2=ssc, 3=depth, 4=lat, 5=long, 6=season
fclose(fid);

yr_ssc_lf=[yr_ssc_lf;T{2}(T{5}>96 & strcmp('Low',T{6}))];
yr_ssc_hf=[yr_ssc_hf;T{2}(T{5}>96 & strcmp('High',T{6}))];

br_ssc_lf=[br_ssc_lf;T{2}(T{5}<96 & T{5}>95 & strcmp('Low',T{6}))];
br_ssc_hf=[br_ssc_hf;T{2}(T{5}<96 & T{5}>95 & strcmp('High',T{6}))];

pr_ssc_lf=[pr_ssc_lf;T{2}(T{5}<95 & strcmp('Low',T{6}))];
pr_ssc_hf=[pr_ssc_hf;T{2}(T{5}<95 & strcmp('High',T{6}))];

%% only compare surface ssc to salinity in casts colocated with water sample
clear all,close all,clc
load('AyeMar18_CTD_all.mat')
cc=1;
for jj=1:length(CombinedProfiles)
    if isempty(CombinedProfiles(jj).SSC)==0
        salinity_surf(cc)=nanmean(CombinedProfiles(jj).Salinity(1:10));
        ssc_surf(cc)=CombinedProfiles(jj).SSC;
        salinity_all(cc)=nanmean(CombinedProfiles(jj).Salinity);
        ssc_all(cc)=nanmean(CombinedProfiles(jj).SSCCal);
        lon(cc)=CombinedProfiles(jj).long;
        cc=cc+1;
    end
end

yr_ssc_lf=[ssc_surf(lon>96),ssc_all(lon>96)];
br_ssc_lf=[ssc_surf(lon<96 & lon>95),ssc_all(lon<96 & lon>95)];
pr_ssc_lf=[ssc_surf(lon<95),ssc_all(lon<95)];

yr_sal_lf=[salinity_surf(lon>96),salinity_all(lon>96)];
br_sal_lf=[salinity_surf(lon<96 & lon>95),salinity_all(lon<96 & lon>95)];
pr_sal_lf=[salinity_surf(lon<95),salinity_all(lon<95)];

load('AyeSept17_CTD_all.mat')
cc=1;
for jj=1:length(CombinedProfiles)
    if isempty(CombinedProfiles(jj).SSC)==0
        salinity_surf(cc)=nanmean(CombinedProfiles(jj).Salinity(1:10));
        ssc_surf(cc)=CombinedProfiles(jj).SSC;
        salinity_all(cc)=nanmean(CombinedProfiles(jj).Salinity);
        ssc_all(cc)=nanmean(CombinedProfiles(jj).SSCCal);
        lon(cc)=CombinedProfiles(jj).long;
        cc=cc+1;
    end
end

yr_ssc_hf=[ssc_surf(lon>96),ssc_all(lon>96)];
br_ssc_hf=[ssc_surf(lon<96 & lon>95),ssc_all(lon<96 & lon>95)];
pr_ssc_hf=[ssc_surf(lon<95),ssc_all(lon<95)];

yr_sal_hf=[salinity_surf(lon>96),salinity_all(lon>96)];
br_sal_hf=[salinity_surf(lon<96 & lon>95),salinity_all(lon<96 & lon>95)];
pr_sal_hf=[salinity_surf(lon<95),salinity_all(lon<95)];

%%
figure;
scatter(yr_ssc_hf,yr_sal_hf,'r*'),hold on
scatter(yr_ssc_lf,yr_sal_lf,'rd')
scatter(br_ssc_hf,br_sal_hf,'b*')
scatter(br_ssc_lf,br_sal_lf,'bd')
scatter(pr_ssc_hf,pr_sal_hf,'k*')
scatter(pr_ssc_lf,pr_sal_lf,'kd')


legend({'Yangon, High','Yangon Low','Bogale High','Bogale Low','Pathein High','Pathein Low'})

ax=gca; ax.XScale='Log';xlim([10 10000])
xlabel('SSC (mg/L)'),ylabel('Salinity')