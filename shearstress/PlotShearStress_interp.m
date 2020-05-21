% plot all of the discharge and bed shear stress figs together

% plot1=discharge and sediment flux
% plot2=shear stress and mean velocity
% plot3 SSC at the thalweg
clear all,close all,clc
cd C:\GLOVER\output\myanmar\
C1=cmocean('turbid');C1=C1(1:200,:);

%%

load('BR_Sept17_Neap_FluxDecomp3.mat')
load('BR_Sept17_Neap_ShearStress.mat')
load('BR_Sept17_Neap_SSC.mat')
% make the ssc matrix
ssc_trans=NaN(length(ssc(1).sigmaAlongComplete(:,1)),length(ssc));
depth=repmat(linspace(0,1,50)',1,length(ssc));
for jj=1:length(ssc)
    depth(:,jj)=depth(:,jj)*max(ssc(jj).adcpbeddepth);
    ssc_trans(:,jj)=nanmean(ssc(jj).sigmaAlongComplete,2);
    ssc_time(jj)=ssc(jj).time(1,1);
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
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Qf,'bo'),hold on
plot(fluxdecomp.time,fluxdecomp.Qf,'b-')
ax=gca;ax.YColor='b';ylim([-11000 11000]),ylabel('Water Flux (m^3/s)')
yyaxis right
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Fs,'k^'),hold on
plot(fluxdecomp.time,fluxdecomp.Fs,'k-')
ax=gca;ax.YColor='k';ylim([-7 7]),ylabel('Sediment Flux (t/s)')
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
colorbar,colormap(C1),caxis([50 350]),ylim([0 20])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])


% Spring
clearvars -except C1

load('BR_Sept17_Spring_FluxDecomp3.mat')
load('BR_Sept17_Spring_ShearStress.mat')
load('BR_Sept17_Spring_SSC.mat')
% make the ssc matrix
ssc_trans=NaN(length(ssc(1).sigmaAlongComplete(:,1)),length(ssc));
depth=repmat(linspace(0,1,50)',1,length(ssc));
for jj=1:length(ssc)
    depth(:,jj)=depth(:,jj)*max(ssc(jj).adcpbeddepth);
    ssc_trans(:,jj)=nanmean(ssc(jj).sigmaAlongComplete,2);
    ssc_time(jj)=ssc(jj).time(1,1);
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
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Qf,'bo'),hold on
plot(fluxdecomp.time,fluxdecomp.Qf,'b-')
ax=gca;ax.YColor='b';ylim([-11000 11000]),ylabel('Water Flux (m^3/s)')
yyaxis right
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Fs,'k^'),hold on
plot(fluxdecomp.time,fluxdecomp.Fs,'k-')
ylabel('Sediment Flux (t/s)')
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
colorbar,colormap(C1),caxis([50 500]),ylim([0 20])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])



%% Bogale River Low Flow:

% plot1 = water level
% plot1=discharge and sediment flux
% plot2=shear stress and mean velocity
% plot3 SSC at the thalweg

clearvars -except C1%clear all,%close all,clc

load('BR_Mar18_Neap_FluxDecomp3.mat')
load('BR_Mar18_Neap_ShearStress.mat')
load('BR_Mar18_Neap_SSC.mat')
% make the ssc matrix
ssc_trans=NaN(length(ssc(1).sigmaAlongComplete(:,1)),length(ssc));
depth=repmat(linspace(0,1,50)',1,length(ssc));
for jj=1:length(ssc)
    depth(:,jj)=depth(:,jj)*max(ssc(jj).adcpbeddepth);
    ssc_trans(:,jj)=nanmean(ssc(jj).sigmaAlongComplete,2);
    ssc_time(jj)=ssc(jj).time(1,1);
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
idx=find(fluxdecomp.time==fluxdecomp.MeasuredTime(end-1));
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Qf,'bo'),hold on
plot(fluxdecomp.time(1:idx),fluxdecomp.Qf(1:idx),'b-')
plot(fluxdecomp.time(idx:end),fluxdecomp.Qf(idx:end),'b--')
ax=gca;ax.YColor='b';ylim([-11000 11000]),ylabel('Water Flux (m^3/s)')
yyaxis right
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Fs,'k^'),hold on
plot(fluxdecomp.time(1:idx),fluxdecomp.Fs(1:idx),'k-')
plot(fluxdecomp.time(idx:end),fluxdecomp.Fs(idx:end),'k--')
ylabel('Sediment Flux (t/s)')
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
colorbar,colormap(C1),caxis([50 500]),ylim([0 20])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])


% Spring
clearvars -except C1

load('BR_Mar18_Spring_FluxDecomp3.mat')
load('BR_Mar18_Spring_ShearStress.mat')
load('BR_Mar18_Spring_SSC.mat')
% make the ssc matrix
ssc_trans=NaN(length(ssc(1).sigmaAlongComplete(:,1)),length(ssc));
depth=repmat(linspace(0,1,50)',1,length(ssc));
for jj=1:length(ssc)
    depth(:,jj)=depth(:,jj)*max(ssc(jj).adcpbeddepth);
    ssc_trans(:,jj)=nanmean(ssc(jj).sigmaAlongComplete,2);
    ssc_time(jj)=ssc(jj).time(1,1);
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
idx=find(fluxdecomp.time==fluxdecomp.MeasuredTime(end-1));
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Qf,'bo'),hold on
plot(fluxdecomp.time(1:idx),fluxdecomp.Qf(1:idx),'b-')
plot(fluxdecomp.time(idx:end),fluxdecomp.Qf(idx:end),'b--')
ax=gca;ax.YColor='b';ylim([-11000 11000]),ylabel('Water Flux (m^3/s)')
yyaxis right
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Fs,'k^'),hold on
plot(fluxdecomp.time(1:idx),fluxdecomp.Fs(1:idx),'k-')
plot(fluxdecomp.time(idx:end),fluxdecomp.Fs(idx:end),'k--')
ylabel('Sediment Flux (t/s)')
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
colorbar,colormap(C1),caxis([50 500]),ylim([0 20])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])





%% Yangon
clearvars -except C1% clear all%,close all,clc

load('YR_Sept17_Neap_FluxDecomp3.mat')
load('YR_Sept17_Neap_ShearStress.mat')
load('YR_Sept17_Neap_SSC.mat')

% make the ssc matrix
ssc_trans=NaN(length(ssc(1).sigmaAlongComplete(:,1)),length(ssc));
depth=repmat(linspace(0,1,50)',1,length(ssc));
for jj=1:length(ssc)
    depth(:,jj)=depth(:,jj)*max(ssc(jj).adcpbeddepth);
    ssc_trans(:,jj)=nanmean(ssc(jj).sigmaAlongComplete,2);
    ssc_time(jj)=ssc(jj).time(1,1);
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
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Qf,'bo'),hold on
plot(fluxdecomp.time,fluxdecomp.Qf,'b-')
ax=gca;ax.YColor='b';ylim([-18000 18000]),ylabel('Water Flux (m^3/s)')
yyaxis right
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Fs,'k^'),hold on
plot(fluxdecomp.time,fluxdecomp.Fs,'k-')
ylabel('Sediment Flux (t/s)')
ax=gca;ax.YColor='k';ylim([-8 8])
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
colorbar,colormap(C1),caxis([50 500]),ylim([0 24])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])


% Yangon neap, Low Flow
clearvars -except C1

load('YR_Jan19_Neap_FluxDecomp3.mat')
load('YR_Jan19_Neap_ShearStress.mat')
load('YR_Jan19_Neap_SSC.mat')

% make the ssc matrix
ssc_trans=NaN(length(ssc(1).sigmaAlongComplete(:,1)),length(ssc));
depth=repmat(linspace(0,1,50)',1,length(ssc));
for jj=1:length(ssc)
    depth(:,jj)=depth(:,jj)*max(ssc(jj).adcpbeddepth);
    ssc_trans(:,jj)=nanmean(ssc(jj).sigmaAlongComplete,2);
    ssc_time(jj)=ssc(jj).time(1,1);
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
idx=find(fluxdecomp.time==fluxdecomp.MeasuredTime(end-1));
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Qf,'bo'),hold on
plot(fluxdecomp.time(1:idx),fluxdecomp.Qf(1:idx),'b-')
plot(fluxdecomp.time(idx:end),fluxdecomp.Qf(idx:end),'b--')
ax=gca;ax.YColor='b';ylim([-18000 18000]),ylabel('Water Flux (m^3/s)')
yyaxis right
plot(fluxdecomp.MeasuredTime,fluxdecomp.Meas.Fs,'k^'),hold on
plot(fluxdecomp.time(1:idx),fluxdecomp.Fs(1:idx),'k-')
plot(fluxdecomp.time(idx:end),fluxdecomp.Fs(idx:end),'k--')
ylabel('Sediment Flux (t/s)')
ax=gca;ax.YColor='k';ylim([-40 40])
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
colorbar,colormap(C1),caxis([0 1500]),ylim([0 24])
xlim(xdates)
datetick('x','HHMM','keeplimits')
ylabel('Depth (m)'),xlabel(['Hour, ',datestr(round(xdates(1)))])
%