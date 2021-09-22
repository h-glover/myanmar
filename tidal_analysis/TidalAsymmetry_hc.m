clear all,close all,clc
cd C:\GLOVER\output\myanmar

load('longterminst\BogaleRiverInstruments.mat'),br=BogaleRiver; clear Bog*

% calculate asymmetry and skew:
% [skew_elev,skew_vel,skew_flux,tvec]=skew_asym(slope,spd,ssc2./1000,time);

% compare 4 day period to 2 week period to full period


%% cyclone High Flow duration asymmetry
hf = 132874:144065;
time = br.datenum(hf)';
% calc the time step in hours:
dt=24*(time(2)-time(1));

days4 = hf(1):round(hf(1)+(4*24)/dt);
weeks2 = hf(1):round(hf(1)+(14*24)/dt);

elev = br.MeinmahlaDepth(hf)-nanmean(br.MeinmahlaDepth(hf));

slope = (elev(3:end) - elev(1:end-2))/(dt*2);
slope = [slope(1),slope,slope(end)];

% moving mean window length for a single tide:
T_move=floor(12.5/dt);

% calculate lenght of time series and make new index vector for sampling at
% 1 : moving mean interval : length of series
L = length(slope);
ind=1:T_move:L; 
tvec = time(ind(1:end-1));

% calculate asymmetry from the water slope: dh/dt
for n = 1:length(ind)-1
    moment_move = slope(ind(n):ind(n+1)-1);
    u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
    sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
    skew_elev_move(n) = u3_move./(sigma3_move);
end

figure;
subplot(411)
yyaxis left,plot(time,elev,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(tvec,skew_elev_move,'k'),hold on
ylabel('\gamma _d'),ylim([-1 1])


nanmean(skew_elev_move)
%%
clearvars -except br hf dt
weeks2 = hf(1):round(hf(1)+(14*24)/dt);

time = br.datenum(weeks2);
elev = br.MeinmahlaDepth(weeks2)-nanmean(br.MeinmahlaDepth(weeks2));

slope = (elev(3:end) - elev(1:end-2))/(dt*2);
slope = [slope(1),slope,slope(end)];

% moving mean window length for a single tide:
T_move=floor(12.5/dt);

% calculate lenght of time series and make new index vector for sampling at
% 1 : moving mean interval : length of series
L = length(slope);
ind=1:T_move:L; 
tvec = time(ind(1:end-1));

% calculate asymmetry from the water slope: dh/dt
for n = 1:length(ind)-1
    moment_move = slope(ind(n):ind(n+1)-1);
    u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
    sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
    skew_elev_move(n) = u3_move./(sigma3_move);
end

yyaxis right,plot(tvec,skew_elev_move,'g--')
xlim([time(1) time(end)])

%%
clearvars -except br hf dt
days4 = hf(1):round(hf(1)+(4*24)/dt);

time = br.datenum(days4);
elev = br.MeinmahlaDepth(days4)-nanmean(br.MeinmahlaDepth(days4));

slope = (elev(3:end) - elev(1:end-2))/(dt*2);
slope = [slope(1),slope,slope(end)];

% moving mean window length for a single tide:
T_move=floor(12.5/dt);

% calculate lenght of time series and make new index vector for sampling at
% 1 : moving mean interval : length of series
L = length(slope);
ind=1:T_move:L; 
tvec = time(ind(1:end-1));

% calculate asymmetry from the water slope: dh/dt
for n = 1:length(ind)-1
    moment_move = slope(ind(n):ind(n+1)-1);
    u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
    sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
    skew_elev_move(n) = u3_move./(sigma3_move);
end

yyaxis right,plot(tvec,skew_elev_move,'r:')
ylabel('\gamma _d'),ylim([-1 1])
refline(0,0)
datetick('x','dd/mm','keeplimits')

legend({'depth','all','2 weeks','4 days'})

%% load aquadopps and compare for the same time periods
clear all,close all,clc
load('aqd\aqd_sep17_hc'); hf17=aqd;
load('aqd\aqd_mar18_hc'); lf18=aqd;
load('aqd\aqd_sep19_hc'); hf19=aqd;clear aqd

% smooth the slope and fill gaps
hf17.slope = fillmissing(hf17.slope,'linear');
hf17.slope = movmean(hf17.slope,5);
hf19.slope = fillmissing(hf19.slope,'linear');
hf19.slope = movmean(hf19.slope,5);
hf19.slope(hf19.slope>1.12)=NaN;

% remove mean from depth for plot
hf17.depth = hf17.depth-nanmean(hf17.depth);
hf19.depth = hf19.depth-nanmean(hf19.depth);
lf18.depth = lf18.depth-nanmean(lf18.depth);

% add sign back to velocity magnitude (ebb negative)
hf17.spd_mean(hf17.slope<0)=(-1)*hf17.spd_mean(hf17.slope<0);
hf19.spd_mean(hf19.slope<0)=(-1)*hf19.spd_mean(hf19.slope<0);
lf18.spd_mean(lf18.slope<0)=(-1)*lf18.spd_mean(lf18.slope<0);

% calculate acceleration
deltaT=round((1/24)./(hf17.time(2)-hf17.time(1)));
hf17.acc=(hf17.spd_mean(3:end)-hf17.spd_mean(1:end-2))./(2*deltaT);
hf17.acc=[hf17.acc(1);hf17.acc;hf17.acc(end)];
hf17.acc(abs(hf17.acc)>0.005)=NaN;
hf17.acc = fillmissing(hf17.acc,'nearest');

% calculate acceleration
deltaT=round((1/24)./(lf18.time(2)-lf18.time(1)));
lf18.acc=(lf18.spd_mean(3:end)-lf18.spd_mean(1:end-2))./(2*deltaT);
lf18.acc=[lf18.acc(1);lf18.acc;lf18.acc(end)];
lf18.acc(abs(lf18.acc)>0.007)=NaN;
lf18.acc = fillmissing(lf18.acc,'nearest');

% calculate acceleration
deltaT=round((1/24)./(hf19.time(2)-hf19.time(1)));
hf19.acc=(hf19.spd_mean(3:end)-hf19.spd_mean(1:end-2))./(2*deltaT);
hf19.acc=[hf19.acc(1);hf19.acc;hf19.acc(end)];
hf19.acc(abs(hf19.acc)>0.0045)=NaN;
hf19.acc = fillmissing(hf19.acc,'nearest');

% calculate asymmetry and skew:
[hf17.skew_elev,hf17.skew_vel,hf17.skew_acc,hf17.tvec]=skew_asym3(hf17.slope,hf17.spd_mean,hf17.acc,hf17.time);
[hf19.skew_elev,hf19.skew_vel,hf19.skew_acc,hf19.tvec]=skew_asym3(hf19.slope,hf19.spd_mean,hf19.acc,hf19.time);
[lf18.skew_elev,lf18.skew_vel,lf18.skew_acc,lf18.tvec]=skew_asym3(lf18.slope,lf18.spd_mean,lf18.acc,lf18.time);

nanmean(hf17.skew_acc)
nanmean(hf19.skew_acc)
nanmean(lf18.skew_acc)

nanmean(hf17.skew_vel)
nanmean(hf19.skew_vel)
nanmean(lf18.skew_vel)

subplot(412),
yyaxis left,plot(hf17.time,hf17.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf17.tvec,hf17.skew_elev,'k')
ylabel('\gamma _d'),ylim([-1 1]),refline(0,0)
xlim([hf17.time(1) hf17.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(413)
yyaxis left,plot(lf18.time,lf18.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(lf18.tvec,lf18.skew_elev,'k')
ylabel('\gamma _d'),ylim([-1 1]),refline(0,0)
xlim([lf18.time(1) lf18.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(414)
yyaxis left,plot(hf19.time,hf19.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf19.tvec,hf19.skew_elev,'k')
ylabel('\gamma _d'),ylim([-1 1]),refline(0,0)
xlim([hf19.time(1) hf19.time(end)]),datetick('x','dd/mm','keeplimits')

figure;
subplot(412),
yyaxis left,plot(hf17.time,hf17.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf17.tvec,hf17.skew_vel,'k')
ylabel('\gamma _v'),ylim([-1 1]),refline(0,0)
xlim([hf17.time(1) hf17.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(413)
yyaxis left,plot(lf18.time,lf18.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(lf18.tvec,lf18.skew_vel,'k')
ylabel('\gamma _v'),ylim([-1 1]),refline(0,0)
xlim([lf18.time(1) lf18.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(414)
yyaxis left,plot(hf19.time,hf19.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf19.tvec,hf19.skew_vel,'k')
ylabel('\gamma _v'),ylim([-1 1]),refline(0,0)
xlim([hf19.time(1) hf19.time(end)]),datetick('x','dd/mm','keeplimits')

figure;
subplot(412),
yyaxis left,plot(hf17.time,hf17.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf17.tvec,hf17.skew_acc,'k')
ylabel('\gamma _v'),ylim([-2 2]),refline(0,0)
xlim([hf17.time(1) hf17.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(413)
yyaxis left,plot(lf18.time,lf18.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(lf18.tvec,lf18.skew_acc,'k')
ylabel('\gamma _v'),ylim([-2 2]),refline(0,0)
xlim([lf18.time(1) lf18.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(414)
yyaxis left,plot(hf19.time,hf19.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf19.tvec,hf19.skew_acc,'k')
ylabel('\gamma _v'),ylim([-2 2]),refline(0,0)
xlim([hf19.time(1) hf19.time(end)]),datetick('x','dd/mm','keeplimits')

figure;
subplot(412),
yyaxis left,plot(hf17.time,hf17.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf17.time,hf17.spd_mean,'k')
ylabel('\gamma _v')%,ylim([-2 2]),refline(0,0)
xlim([hf17.time(1) hf17.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(413)
yyaxis left,plot(lf18.time,lf18.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(lf18.time,lf18.spd_mean,'k')
ylabel('\gamma _v'),ylim([-2 2]),refline(0,0)
xlim([lf18.time(1) lf18.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(414)
yyaxis left,plot(hf19.time,hf19.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf19.time,hf19.spd_mean,'k')
ylabel('\gamma _v'),ylim([-2 2]),refline(0,0)
xlim([hf19.time(1) hf19.time(end)]),datetick('x','dd/mm','keeplimits')


