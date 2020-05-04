% plot all of the discharge and bed shear stress figs together

% plot1 = water level
% plot1=discharge and sediment flux
% plot2=shear stress and mean velocity
% plot3 SSC at the thalweg

clear all,close all,clc

load('BR_Sept17_Neap_FluxDecomp2.mat')
load('BR_Sept17_Neap_ShearStress.mat')
load('BR_Sept17_Neap_SSC.mat')
% make the ssc matrix
depth=0:0.1:18;
ssc_trans=NaN(length(depth),length(ssc));
for jj=1:length(ssc)
    [dd,~]=max(ssc(jj).Depth);
    [~,idx]=max(dd); 
    ssc(jj).Depth(:,idx)=fillmissing(ssc(jj).Depth(:,idx),'linear');
    ssc_trans(:,jj)=interp1(ssc(jj).Depth(:,idx),ssc(jj).SSC(:,idx),depth);
    ssc_time(jj)=ssc(jj).time(1,idx);
end

% take the avg of all the tau
tau_qsl=cat(3,avgs.tau_qsl);
tau_qsl_std=std(tau_qsl,0,3);
tau_qsl=nanmean(tau_qsl,3);
%tau_time=repmat(avgs(1).time',1,3);
tau_time=avgs(1).time;

% set same date range for all plots
xdates=[fluxdecomp.time(1)-(1/48) fluxdecomp.time(end)+(1/48)];

% make flood velocities negative
negs=fluxdecomp.time(fluxdecomp.Qf<0);
[~,inegs,itau]=intersect(negs,tau_time);
ubar=avgs(1).ubar(:,2)./100;
ubar(itau)=(-1)*ubar(itau);

figure;
subplot(321)% discharge
yyaxis left
plot(fluxdecomp.time,fluxdecomp.Qf,'b'),ylabel('Water Flux (m^3/s)')
ax=gca;ax.YColor='b';ylim([-11000 11000])
yyaxis right
plot(fluxdecomp.time,fluxdecomp.Fs,'k'),ylabel('Sediment Flux (t/s)')
ax=gca;ax.YColor='k';ylim([-7 7])
xlim(xdates)
datetick('x','HHMM','keeplimits')
title('Bogale River, High Flow, Neap')

subplot(323)% bed stress
yyaxis left
errorbar(tau_time,tau_qsl(:,2),tau_qsl_std(:,2),'r'),ylabel('Tau (Pa)')
ax=gca;ax.YColor='r';ylim([0 7])
yyaxis right
plot(tau_time,ubar,'k'),ylabel('Mean Velocity (m/s)')
ax=gca;ax.YColor='k';ylim([-1.1 1.1])
xlim(xdates)
datetick('x','HHMM','keeplimits')

% ssc
subplot(325)
pcolor(ssc_time,depth,ssc_trans),shading interp,axis ij
colorbar,caxis([0 500]),ylim([0 16])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])


% Spring
clear all

load('BR_Sept17_Spring_FluxDecomp2.mat')
load('BR_Sept17_Spring_ShearStress.mat')
load('BR_Sept17_Spring_SSC.mat')
% make the ssc matrix
depth=0:0.1:18;
ssc_trans=NaN(length(depth),length(ssc));
for jj=1:length(ssc)
    [dd,~]=max(ssc(jj).Depth);
    [~,idx]=max(dd); 
    ssc(jj).Depth(:,idx)=fillmissing(ssc(jj).Depth(:,idx),'linear');
    ssc_trans(:,jj)=interp1(ssc(jj).Depth(:,idx),ssc(jj).SSC(:,idx),depth);
    ssc_time(jj)=ssc(jj).time(1,idx);
end

% take the avg of all the tau
tau_qsl=cat(3,avgs.tau_qsl);
tau_qsl_std=std(tau_qsl,0,3);
tau_qsl=nanmean(tau_qsl,3);
%tau_time=repmat(avgs(1).time',1,3);
tau_time=datevec(avgs(1).time);tau_time(:,6)=0;tau_time=datenum(tau_time);

% set same date range for all plots
xdates=[fluxdecomp.time(1)-(1/48) fluxdecomp.time(end)+(1/48)];

% make flood velocities negative
negs=datevec(fluxdecomp.time(fluxdecomp.Qf<0));negs(:,6)=0;negs=datenum(negs);
[~,~,itau]=intersect(negs,tau_time);
ubar=avgs(1).ubar(:,2)./100;
ubar(itau)=(-1)*ubar(itau);

%figure;
subplot(322)% discharge
yyaxis left
plot(fluxdecomp.time,fluxdecomp.Qf,'b'),ylabel('Water Flux (m^3/s)')
ax=gca;ax.YColor='b';ylim([-11000 11000])
yyaxis right
plot(fluxdecomp.time,fluxdecomp.Fs,'k'),ylabel('Sediment Flux (t/s)')
ax=gca;ax.YColor='k';ylim([-7 7])
xlim(xdates)
datetick('x','HHMM','keeplimits')
title('Bogale River, High Flow, Spring')

subplot(324)% bed stress
yyaxis left
errorbar(tau_time,tau_qsl(:,2),tau_qsl_std(:,2),'r'),ylabel('Tau (Pa)')
ax=gca;ax.YColor='r';ylim([0 7])
yyaxis right
plot(tau_time,ubar,'k'),ylabel('Mean Velocity (m/s)')
ax=gca;ax.YColor='k';ylim([-1.1 1.1])
xlim(xdates)
datetick('x','HHMM','keeplimits')

% ssc
subplot(326)
pcolor(ssc_time,depth,ssc_trans),shading interp,axis ij
colorbar,caxis([0 500]),ylim([0 16])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])



% Bogale River Low Flow:

% plot1 = water level
% plot1=discharge and sediment flux
% plot2=shear stress and mean velocity
% plot3 SSC at the thalweg

clear all%,close all,clc

load('BR_Mar18_Neap_FluxDecomp2.mat')
load('BR_Mar18_Neap_ShearStress.mat')
load('BR_Mar18_Neap_SSC.mat')
% make the ssc matrix
depth=0:0.1:18;
ssc_trans=NaN(length(depth),length(ssc));
for jj=1:length(ssc)
    [dd,~]=max(ssc(jj).Depth);
    [~,idx]=max(dd); 
    ssc(jj).Depth(:,idx)=fillmissing(ssc(jj).Depth(:,idx),'linear');
    ssc_trans(:,jj)=interp1(ssc(jj).Depth(:,idx),ssc(jj).SSC(:,idx),depth);
    ssc_time(jj)=ssc(jj).time(1,idx);
end

% take the avg of all the tau
tau_qsl=cat(3,avgs.tau_qsl);
tau_qsl_std=std(tau_qsl,0,3,'omitnan');
tau_qsl=nanmean(tau_qsl,3);
%tau_time=repmat(avgs(1).time',1,3);
tau_time=datevec(avgs(1).time);tau_time(:,6)=0;tau_time=datenum(tau_time);

% set same date range for all plots
xdates=[fluxdecomp.time(1)-(1/48) fluxdecomp.time(end)+(1/48)];

% make flood velocities negative
negs=datevec(fluxdecomp.time(fluxdecomp.Qf<0));negs(:,6)=0;negs=datenum(negs);
[~,~,itau]=intersect(negs,tau_time);
ubar=avgs(1).ubar(:,2)./100;
ubar(itau)=(-1)*ubar(itau);


figure;
subplot(321)% discharge
yyaxis left
plot(fluxdecomp.time,fluxdecomp.Qf,'b'),ylabel('Water Flux (m^3/s)')
ax=gca;ax.YColor='b';ylim([-11000 11000])
yyaxis right
plot(fluxdecomp.time,fluxdecomp.Fs,'k'),ylabel('Sediment Flux (t/s)')
ax=gca;ax.YColor='k';ylim([-7 7])
xlim(xdates)
datetick('x','HHMM','keeplimits')
title('Bogale River, Low Flow, Neap')

subplot(323)% bed stress
yyaxis left
errorbar(tau_time,tau_qsl(:,2),tau_qsl_std(:,2),'r'),ylabel('Tau (Pa)')
ax=gca;ax.YColor='r';ylim([0 7])
yyaxis right
plot(tau_time,ubar,'k'),ylabel('Mean Velocity (m/s)')
ax=gca;ax.YColor='k';ylim([-1.1 1.1])
xlim(xdates)
datetick('x','HHMM','keeplimits')

% ssc
subplot(325)
pcolor(ssc_time,depth,ssc_trans),shading interp,axis ij
colorbar,caxis([0 500]),ylim([0 16])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])


% Spring
clear all

load('BR_Mar18_Spring_FluxDecomp2.mat')
load('BR_Mar18_Spring_ShearStress.mat')
load('BR_Mar18_Spring_SSC.mat')
% make the ssc matrix
depth=0:0.1:18;
ssc_trans=NaN(length(depth),length(ssc));
for jj=1:length(ssc)
    [dd,~]=max(ssc(jj).Depth);
    [~,idx]=max(dd); 
    ssc(jj).Depth(:,idx)=fillmissing(ssc(jj).Depth(:,idx),'linear');
    ssc_trans(:,jj)=interp1(ssc(jj).Depth(:,idx),ssc(jj).SSC(:,idx),depth);
    ssc_time(jj)=ssc(jj).time(1,idx);
end

% take the avg of all the tau
tau_qsl=cat(3,avgs.tau_qsl);
tau_qsl_std=std(tau_qsl,0,3,'omitnan');
tau_qsl=nanmean(tau_qsl,3);
%tau_time=repmat(avgs(1).time',1,3);
tau_time=datevec(avgs(1).time);tau_time(:,6)=0;tau_time=datenum(tau_time);

% set same date range for all plots
xdates=[fluxdecomp.time(1)-(1/48) fluxdecomp.time(end)+(1/48)];

% make flood velocities negative
negs=datevec(fluxdecomp.time(fluxdecomp.Qf<0));negs(:,6)=0;negs=datenum(negs);
[~,~,itau]=intersect(negs,tau_time);
ubar=avgs(1).ubar(:,2)./100;
ubar(itau)=(-1)*ubar(itau);




%figure;
subplot(322)% discharge
yyaxis left
plot(fluxdecomp.time,fluxdecomp.Qf,'b'),ylabel('Water Flux (m^3/s)')
ax=gca;ax.YColor='b';ylim([-11000 11000])
yyaxis right
plot(fluxdecomp.time,fluxdecomp.Fs,'k'),ylabel('Sediment Flux (t/s)')
ax=gca;ax.YColor='k';ylim([-7 7])
xlim(xdates)
datetick('x','HHMM','keeplimits')
title('Bogale River, Low Flow, Spring')

subplot(324)% bed stress
yyaxis left
errorbar(tau_time,tau_qsl(:,2),tau_qsl_std(:,2),'r'),ylabel('Tau (Pa)')
ax=gca;ax.YColor='r';ylim([0 7])
yyaxis right
plot(tau_time,ubar,'k'),ylabel('Mean Velocity (m/s)')
ax=gca;ax.YColor='k';ylim([-1.1 1.1])
xlim(xdates)
datetick('x','HHMM','keeplimits')

% ssc
subplot(326)
pcolor(ssc_time,depth,ssc_trans),shading interp,axis ij
colorbar,caxis([0 500]),ylim([0 16])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])





%% Yangon
clear all,close all,clc

load('YR_Sept17_Neap_FluxDecomp2.mat')
load('YR_Sept17_Neap_ShearStress.mat')
load('YR_Sept17_Neap_SSC.mat')

% make the ssc matrix
depth=0:0.1:24;
ssc_trans=NaN(length(depth),length(ssc));
for jj=1:length(ssc)
    [dd,~]=max(ssc(jj).Depth);
    [~,idx]=max(dd); 
    ssc(jj).Depth(:,idx)=fillmissing(ssc(jj).Depth(:,idx),'linear');
    ssc_trans(:,jj)=interp1(ssc(jj).Depth(:,idx),ssc(jj).SSC(:,idx),depth);
    ssc_time(jj)=ssc(jj).time(1,idx);
end

% take the avg of all the tau
tau_qsl=cat(3,avgs.tau_qsl);
tau_qsl_std=std(tau_qsl,0,3);
tau_qsl=nanmean(tau_qsl,3);
%tau_time=repmat(avgs(1).time',1,3);
tau_time=datevec(avgs(1).time);tau_time(:,6)=0;tau_time=datenum(tau_time);

% set same date range for all plots
xdates=[fluxdecomp.time(1)-(1/48) fluxdecomp.time(end)+(1/48)];

% make flood velocities negative
negs=datevec(fluxdecomp.time(fluxdecomp.Qf<0));negs(:,6)=0;negs=datenum(negs);
[~,~,itau]=intersect(negs,tau_time);
ubar=avgs(1).ubar(:,2)./100;
ubar(itau)=(-1)*ubar(itau);


figure;
subplot(321)% discharge
yyaxis left
plot(fluxdecomp.time,fluxdecomp.Qf,'b'),ylabel('Water Flux (m^3/s)')
ax=gca;ax.YColor='b';ylim([-11000 18000])
yyaxis right
plot(fluxdecomp.time,fluxdecomp.Fs,'k'),ylabel('Sediment Flux (t/s)')
ax=gca;ax.YColor='k';ylim([-2 8])
xlim(xdates)
datetick('x','HHMM','keeplimits')
title('Yangon River, High Flow, Neap')

subplot(323)% bed stress
yyaxis left
errorbar(tau_time,tau_qsl(:,2),tau_qsl_std(:,2),'r'),ylabel('Tau (Pa)')
ax=gca;ax.YColor='r';ylim([0 12])
yyaxis right
plot(tau_time,ubar,'k'),ylabel('Mean Velocity (m/s)')
ax=gca;ax.YColor='k';ylim([-1.1 1.6])
xlim(xdates)
datetick('x','HHMM','keeplimits')

% ssc
subplot(325)
pcolor(ssc_time,depth,ssc_trans),shading interp,axis ij
colorbar,caxis([0 600]),ylim([0 18])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])


% Yangon neap, Low Flow
clear all%,close all,clc

load('YR_Jan19_Neap_FluxDecomp2.mat')
load('YR_Jan19_Neap_ShearStress.mat')
load('YR_Jan19_Neap_SSC.mat')

% make the ssc matrix
depth=0:0.1:24;
ssc_trans=NaN(length(depth),length(ssc));
for jj=1:length(ssc)
    [dd,~]=max(ssc(jj).Depth);
    [~,idx]=max(dd); 
    ssc(jj).Depth(:,idx)=fillmissing(ssc(jj).Depth(:,idx),'linear');
    ssc_trans(:,jj)=interp1(ssc(jj).Depth(:,idx),ssc(jj).SSC(:,idx),depth);
    ssc_time(jj)=ssc(jj).time(1,idx);
end

% take the avg of all the tau
tau_qsl=cat(3,avgs.tau_qsl);
tau_qsl_std=std(tau_qsl,0,3,'omitnan');
tau_qsl=nanmean(tau_qsl,3);
%tau_time=repmat(avgs(1).time',1,3);
tau_time=datevec(avgs(1).time);tau_time(:,6)=0;tau_time=datenum(tau_time);

% set same date range for all plots
xdates=[fluxdecomp.time(1)-(1/48) fluxdecomp.time(end)+(1/48)];

% make flood velocities negative
negs=datevec(fluxdecomp.time(fluxdecomp.Qf<0));negs(:,6)=0;negs=datenum(negs);
[~,~,itau]=intersect(negs,tau_time);
ubar=avgs(1).ubar(:,3)./100;
ubar(itau)=(-1)*ubar(itau);

subplot(322)% discharge
yyaxis left
plot(fluxdecomp.time,fluxdecomp.Qf,'b'),ylabel('Water Flux (m^3/s)')
ax=gca;ax.YColor='b';ylim([-11000 18000])
yyaxis right
plot(fluxdecomp.time,fluxdecomp.Fs,'k'),ylabel('Sediment Flux (t/s)')
ax=gca;ax.YColor='k';ylim([-28 34])
xlim(xdates)
datetick('x','HHMM','keeplimits')
title('Yangon River, Low Flow, Neap')

subplot(324)% bed stress
yyaxis left
errorbar(tau_time,tau_qsl(:,3),tau_qsl_std(:,3),'r'),ylabel('Tau (Pa)')
ax=gca;ax.YColor='r';ylim([0 12])
yyaxis right
plot(tau_time,ubar,'k'),ylabel('Mean Velocity (m/s)')
ax=gca;ax.YColor='k';ylim([-1.1 1.6])
xlim(xdates)
datetick('x','HHMM','keeplimits')

% ssc
subplot(326)
pcolor(ssc_time,depth,ssc_trans),shading interp,axis ij
colorbar,caxis([0 1500]),ylim([0 18])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])
 %% 
% % create a height-above-bed vector for the new bin coordinates (cm)
% hab_mat=100:10:1800;
% hab_log=log(hab_mat);
% zdepth=round(adcp(1).z(:,1).*100);
% xbins=20;
% for jj=1:length(adcp)
%     % calculate the cross-section oriented speed(cm/s):
%     adcp(jj).interpSpeed=sqrt(adcp(jj).interpalong.^2 + adcp(jj).interpacross.^2);
%     
%     % fix interp depth and convert to cm
%     adcp(jj).interpdepths(...
%         adcp(jj).interpdepths<0 | adcp(jj).interpdepths>26)=NaN;
%     adcp(jj).interpdepths=round(adcp(jj).interpdepths.*100);
% 
%     % calculate the actual height above bed of each profile's bins (cm)
%     adcp(jj).hab_actual=adcp(jj).interpdepths-zdepth;
%     
%     % for each column with depth>0, interpolate to the new HAB matrix so
%     % that all of the profile bins line up
%     cols=find(~isnan(adcp(jj).interpdepths));
%     adcp(jj).spd_hab=NaN([length(hab_mat),length(adcp(jj).interpdepths)]);
%     for nn=cols
%         adcp(jj).spd_hab(:,nn)=interp1(...
%             adcp(jj).hab_actual(:,nn),adcp(jj).interpSpeed(:,nn),hab_mat);
%     end
%     MeasuredTime(jj)=adcp(jj).time(1);
%     spd_hab(:,jj)=nanmedian(adcp(jj).spd_hab(:,300-xbins:300+xbins),2);
% end