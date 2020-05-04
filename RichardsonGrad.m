% calculating Ri for changes in density
clear all,close all,clc
load('AyeMar18_CTD_all.mat')
load('BR_Mar18_Spring.mat')
ctd=CombinedProfiles;clear CombinedProfiles
for jj=1:length(ctd)
    time(jj)=ctd(jj).time(1);
end

time=datevec(time); 
ctd=ctd(time(:,3)==5);

%
nn=5;
cidx=repmat([1/nn],1,nn);

depth=0:0.1:20;


for jj=1:length(ctd)
    ctd(jj).Density_smooth=conv(ctd(jj).Density,cidx,'same');
    ctd(jj).Density_smooth([1:3,end-3:end])=NaN;
    ctd(jj).Depth_smooth=conv(ctd(jj).Depth,cidx,'same');
    ctd(jj).Depth_smooth([1:nn,end-nn:end])=ctd(jj).Depth([1:nn,end-nn:end]);
    ctd(jj).Depth_smooth=fillmissing(ctd(jj).Depth_smooth,'linear');
    
    ctd(jj).Density_interp=interp1(...
        ctd(jj).Depth_smooth,ctd(jj).Density_smooth,depth);
    ctd(jj).deltarho=(ctd(jj).Density_interp(2:end)-ctd(jj).Density_interp(1:end-1))...
        /0.1;
    figure;
    subplot(121)
    plot(ctd(jj).Density,ctd(jj).Depth,'k'),hold on
    plot(ctd(jj).Density_smooth,ctd(jj).Depth_smooth,'r')
    plot(ctd(jj).Density_interp,depth,'o'),axis ij
    subplot(122)
    plot(ctd(jj).deltarho,depth(2:end)),axis ij,hold on
end
