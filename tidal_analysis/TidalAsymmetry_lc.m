clear all,close all,clc
cd C:\GLOVER\output\myanmar

load('aqd\aqd_sep17_lc'); hf17=aqd;
load('aqd\aqd_mar18_lc'); lf18=aqd;

% figure;
% subplot(211)
% yyaxis left,plot(hf17.depth)
% yyaxis right,plot(hf17.spd_mean)
% subplot(212)
% yyaxis left,plot(hf17.depth)
% yyaxis right,plot(hf17.ssc1)
% 
% figure;
% subplot(211)
% yyaxis left,plot(lf18.depth)
% yyaxis right,plot(lf18.spd_mean)
% subplot(212)
% yyaxis left,plot(lf18.depth)
% yyaxis right,plot(lf18.ssc1)

%%
% smooth the slope and fill gaps
hf17.slope = fillmissing(hf17.slope,'linear');
hf17.slope = movmean(hf17.slope,5);
hf17.spd_mean = movmean(hf17.spd_mean,5);
lf18.spd_mean = movmean(lf18.spd_mean,5);

% remove mean from depth for plot
hf17.depth = hf17.depth-nanmean(hf17.depth);
lf18.depth = lf18.depth-nanmean(lf18.depth);

% add sign back to velocity magnitude (ebb negative)
hf17.spd_mean(hf17.slope<0)=(-1)*hf17.spd_mean(hf17.slope<0);
lf18.spd_mean(lf18.slope<0)=(-1)*lf18.spd_mean(lf18.slope<0);

% calculate acceleration
deltaT=round((1/24)./(hf17.time(2)-hf17.time(1)));
hf17.acc=(hf17.spd_mean(3:end)-hf17.spd_mean(1:end-2))./(2*deltaT);
hf17.acc=[hf17.acc(1);hf17.acc;hf17.acc(end)];
hf17.acc(abs(hf17.acc)>0.0007)=NaN;
hf17.acc = fillmissing(hf17.acc,'nearest');

% calculate acceleration
deltaT=round((1/24)./(lf18.time(2)-lf18.time(1)));
lf18.acc=(lf18.spd_mean(3:end)-lf18.spd_mean(1:end-2))./(2*deltaT);
lf18.acc=[lf18.acc(1);lf18.acc;lf18.acc(end)];
lf18.acc(abs(lf18.acc)>0.015)=NaN;
lf18.acc = fillmissing(lf18.acc,'nearest');

%% calculate asymmetry and skew:
[hf17.skew_elev,hf17.skew_vel,hf17.skew_acc,hf17.tvec]=skew_asym3(hf17.slope,hf17.spd_mean,hf17.acc,hf17.time);
[lf18.skew_elev,lf18.skew_vel,lf18.skew_acc,lf18.tvec]=skew_asym3(lf18.slope,lf18.spd_mean,lf18.acc,lf18.time);

figure;
subplot(211),
yyaxis left,plot(hf17.time,hf17.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf17.tvec,hf17.skew_elev,'k')
ylabel('\gamma _d'),ylim([-1 1]),refline(0,0)
xlim([hf17.time(1) hf17.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(212)
yyaxis left,plot(lf18.time,lf18.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(lf18.tvec,lf18.skew_elev,'k')
ylabel('\gamma _d'),ylim([-1 1]),refline(0,0)
xlim([lf18.time(1) lf18.time(end)]),datetick('x','dd/mm','keeplimits')

figure;
subplot(211),
yyaxis left,plot(hf17.time,hf17.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf17.tvec,hf17.skew_vel,'k')
ylabel('\gamma _v'),ylim([-1 1]),refline(0,0)
xlim([hf17.time(1) hf17.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(212)
yyaxis left,plot(lf18.time,lf18.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(lf18.tvec,lf18.skew_vel,'k')
ylabel('\gamma _v'),ylim([-1 1]),refline(0,0)
xlim([lf18.time(1) lf18.time(end)]),datetick('x','dd/mm','keeplimits')

figure;
subplot(211),
yyaxis left,plot(hf17.time,hf17.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf17.tvec,hf17.skew_acc,'k')
ylabel('\gamma _F'),ylim([-2 2]),refline(0,0)
xlim([hf17.time(1) hf17.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(212)
yyaxis left,plot(lf18.time,lf18.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(lf18.tvec,lf18.skew_acc,'k')
ylabel('\gamma _F'),ylim([-2 2]),refline(0,0)
xlim([lf18.time(1) lf18.time(end)]),datetick('x','dd/mm','keeplimits')

figure;
subplot(211),
yyaxis left,plot(hf17.time,hf17.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(hf17.time,hf17.spd_mean,'k')
ylabel('spd'),ylim([-2 2]),refline(0,0)
xlim([hf17.time(1) hf17.time(end)]),datetick('x','dd/mm','keeplimits')
subplot(212)
yyaxis left,plot(lf18.time,lf18.depth,'b'),hold on
ylabel('depth'),ylim([-2 2])
yyaxis right,plot(lf18.time,lf18.spd_mean,'k')
ylabel('spd'),ylim([-2 2]),refline(0,0)
xlim([lf18.time(1) lf18.time(end)]),datetick('x','dd/mm','keeplimits')


