% compare skew and asym along estuary from observational data:
clear all,close all,clc
cd C:\GLOVER\output\myanmar

load('longterminst\BogaleRiverInstruments.mat'),br=BogaleRiver; clear Bog*

% calculate asymmetry and skew:
% [skew_elev,skew_vel,skew_flux,tvec]=skew_asym(slope,spd,ssc2./1000,time);

%% cyclone High Flow duration asymmetry
hf = 132874:144065;
time = br.datenum(hf)';
elev = br.MeinmahlaDepth(hf)-nanmean(br.MeinmahlaDepth(hf));

% coef = ut_solv (time,elev,[],16,'auto'); 

% calc the time step in hours:
dt=24*(time(2)-time(1));

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
% % Add last point with backward looking window
% moment_move = slope(L-T_move:L);
% u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
% sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
% skew_elev_move(end+1) = u3_move./(sigma3_move);

figure;
subplot(221)
yyaxis left,plot(time,elev,'b')
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(tvec,skew_elev_move,'k')
ylabel('\gamma _d'),ylim([-1 1])
refline(0,0)
title('cyclone')
xlim([br.datenum(hf(1)) br.datenum(hf(end))])
datetick('x','mm','keeplimits')
xlabel('hf, mm')
nanmean(skew_elev_move)
%% cyclone low Flow duration asymmetry
clearvars -except br

lf = 98000:108778;
time = br.datenum(lf)';
elev = br.MeinmahlaDepth(lf)-nanmean(br.MeinmahlaDepth(lf));

% coef = ut_solv (time,elev,[],16,'auto'); 

% calc the time step in hours:
dt=24*(time(2)-time(1));

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

% figure;
subplot(222)
yyaxis left,plot(time,elev,'b')
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(tvec,skew_elev_move,'k')
ylabel('\gamma _d'),ylim([-1 1])
refline(0,0)
title('cyclone')
xlim([br.datenum(lf(1)) br.datenum(lf(end))])
datetick('x','mm','keeplimits')
xlabel('lf, mm')
nanmean(skew_elev_move)
%% FREDA high flow duration asymmetry:
clearvars -except br

% calculate asymmetry and skew:
% [skew_elev,skew_vel,skew_flux,tvec]=skew_asym(slope,spd,ssc2./1000,time);
hf = 133388:143953;
time = br.datenum(hf)';
% elev = br.FredaDepth;
% elev(isnan(elev)) = br.ut_FredaDepth(isnan(elev));
% elev = elev(hf)-nanmean(elev(hf));
elev = br.ut_FredaDepth(hf)-nanmean(br.ut_FredaDepth(hf));

% coef = ut_solv (time,elev,[],16,'auto'); 

% calc the time step in hours:
dt=24*(time(2)-time(1));

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
% % Add last point with backward looking window
% moment_move = slope(L-T_move:L);
% u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
% sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
% skew_elev_move(end+1) = u3_move./(sigma3_move);

subplot(223)
yyaxis left,plot(time,elev,'b')
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(tvec,skew_elev_move,'k')
ylabel('\gamma _d'),ylim([-1 1])
refline(0,0)
title('freda')
xlim([br.datenum(hf(1)) br.datenum(hf(end))])
datetick('x','mm','keeplimits')
xlabel('hf, mm')
nanmean(skew_elev_move)
%% Freda Low Flow duration asymmetry:
clearvars -except br

lf = 45370:55406;
time = br.datenum(lf)';
elev = br.ut_MeinmahlaDepth(lf)-nanmean(br.ut_MeinmahlaDepth(lf));

% coef = ut_solv (time,elev,[],16,'auto'); 

% calc the time step in hours:
dt=24*(time(2)-time(1));

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

% figure;
subplot(224)
yyaxis left,plot(time,elev,'b')
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(tvec,skew_elev_move,'k')
ylabel('\gamma _d'),ylim([-1 1])
refline(0,0)
title('freda')
xlim([br.datenum(lf(1)) br.datenum(lf(end))])
datetick('x','mm','keeplimits')
xlabel('lf, mm')

nanmean(skew_elev_move)
%%  Calculate velocity skewness 
% for n = 1:length(ind)-1
%     moment_move = vel(ind(n):ind(n+1)-1);
%     u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
%     sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
%     skew_vel_move(n) = u3_move./(sigma3_move);
% end
% % Add last point with backward looking window
% moment_move = vel(L-T_move:L);
% u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
% sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
% skew_vel_move(end+1) = u3_move./(sigma3_move);

