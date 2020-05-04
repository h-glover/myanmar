clear all
load YR_Jan19_Neap.mat
for ii=1:length(adcp)
MeasuredDis(ii)=nansum(nansum(adcp(ii).alongComplete./100)).*0.5;

MeasuredTime(ii)=mean(adcp(ii).time);
end
% old, linear interpolation code
% InterpTimes_linear=days([min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.4)]);
% InterpDis_linear=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],...
%     [MeasuredDis,MeasuredDis(1)],InterpTimes_linear,'linear');
% InterpRes=mean(InterpDis_linear);
InterpTimes_linear=days([min(MeasuredTime):minutes(1):max(MeasuredTime)]);
InterpDis_linear=interp1(MeasuredTime,MeasuredDis,InterpTimes_linear,'linear');

% the fill gaps will have interptimes nans interptimes
filltimes=days([MeasuredTime(1):minutes(1):MeasuredTime(1)+hours(12.4)]);
InterpTimes_fg=[filltimes,days([filltimes(end)+minutes(1):minutes(1):...
    filltimes(end)+minutes(1)+hours(12.4)])];
L=length(filltimes)-length(InterpTimes_linear);
InterpDis_fg=[InterpDis_linear,NaN([1,L]),InterpDis_linear,NaN([1,L])];

InterpDis_fg=fillmissing(InterpDis_fg,'spline');
InterpTimes=InterpTimes_fg(1:end-length(filltimes));
InterpDis=InterpDis_fg(1:end-length(filltimes));
InterpRes=mean(InterpDis);

clear adcp


save('YR_neap_discharge_interp','InterpDis','InterpTimes','InterpRes')
%%
clear all
load YR_Jan19_Neap_SedFlux.mat
for ii=1:length(adcp)
MeasuredTime(ii)=mean(adcp(ii).time);
end
MeasuredDis=horzcat(adcp.MeasuredDis);
% InterpTimes=days([min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.4)]);
% InterpDis=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],[MeasuredDis,MeasuredDis(1)],InterpTimes,'linear');
% InterpRes=mean(InterpDis);
% 
MeasuredFlux=horzcat(adcp.sigmaMeasuredSedFlux);
% InterpFlux=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],[MeasuredFlux,MeasuredFlux(1)],InterpTimes,'linear');
% InterpResFlux=mean(InterpFlux);
% 
InterpTimes_linear=days([min(MeasuredTime):minutes(1):max(MeasuredTime)]);
InterpDis_linear=interp1(MeasuredTime,MeasuredDis,InterpTimes_linear,'linear');
InterpFlux_linear=interp1(MeasuredTime,MeasuredFlux,InterpTimes_linear,'linear');

% the fill gaps will have interptimes nans interptimes
filltimes=days([MeasuredTime(1):minutes(1):MeasuredTime(1)+hours(12.4)]);
InterpTimes_fg=[filltimes,days([filltimes(end)+minutes(1):minutes(1):...
    filltimes(end)+minutes(1)+hours(12.4)])];
L=length(filltimes)-length(InterpTimes_linear);
InterpDis_fg=[InterpDis_linear,NaN([1,L]),InterpDis_linear,NaN([1,L])];
InterpFlux_fg=[InterpFlux_linear,NaN([1,L]),InterpFlux_linear,NaN([1,L])];


InterpDis_fg=fillmissing(InterpDis_fg,'spline');
InterpFlux_fg=fillmissing(InterpFlux_fg,'spline');
InterpTimes=InterpTimes_fg(1:end-length(filltimes));
InterpDis=InterpDis_fg(1:end-length(filltimes));
InterpFlux=InterpFlux_fg(1:end-length(filltimes));

InterpRes=mean(InterpDis);
InterpResFlux=mean(InterpFlux);

clear adcp


save('YR_Jan19_Neap_discharge_interp','Interp*','Measured*')
%%
load('YangonWL_Prediction.mat')

figure;
scatter(MeasuredTime,MeasuredFlux,50,'b')
hold on
plot(InterpTimes,InterpFlux,':b','linewidth',2)
set(gca,'fontsize',16)
xlabel('Time','fontsize',16)
ylabel('Discharge (m^3 s^{-1})','fontsize',16,'color','b') % left y-axis 
%ylim([-15000 15000])
%plot(InterpTimes,ones(length(InterpTimes)).*InterpRes,'--b','linewidth',2)

%t=text(mean(InterpTimes),InterpRes-1500,strcat(num2str(round(InterpRes)),' m^3 s^-^1'));
t.FontSize=16;
t.Color='b';

yyaxis right
plot(p_time,YOUT2,'r')
ylim([-3 3])
r=refline(0,0);r.Color='r';
ylabel('Water level')
xlim([min(InterpTimes)-1/48 max(InterpTimes)+1/48])
datetick('x','keeplimits')
