load('AyeMar18_CTD_all.mat')

for jj=1:length(CombinedProfiles)
    CombinedProfiles(jj).Density=sw_dens(CombinedProfiles(jj).Salinity,...
        CombinedProfiles(jj).Temperature,CombinedProfiles(jj).Depth);
    CombinedProfiles(jj).Density_avg=nanmean(CombinedProfiles(jj).Density);

end
subplot(131)
plot(vertcat(CombinedProfiles.Density_avg))
save('AyeMar18_CTD_all.mat','CombinedProfiles')
clear all

load('AyeSept17_CTD_all.mat')
for jj=1:length(CombinedProfiles)
    CombinedProfiles(jj).Density=sw_dens(CombinedProfiles(jj).Salinity,...
        CombinedProfiles(jj).Temperature,CombinedProfiles(jj).Depth);
    CombinedProfiles(jj).Density_avg=nanmean(CombinedProfiles(jj).Density);

end
subplot(132)
plot(vertcat(CombinedProfiles.Density_avg))
save('AyeSept17_CTD_all.mat','CombinedProfiles')
clear all

load('AyeJan19_CTD_all.mat')
for jj=1:length(profiles)
    profiles(jj).Density=sw_dens(profiles(jj).Salinity,...
        profiles(jj).Temperature,profiles(jj).Depth);
    profiles(jj).Density_avg=nanmean(profiles(jj).Density);

end
subplot(133)
plot(vertcat(profiles.Density_avg))
save('AyeJan19_CTD_all.mat','profiles')
