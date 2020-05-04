%% calculate total discharge
clear all
load('') %load adcp data
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

% Use a non-linear fill if you're missing more than ~30mins, otherwise just
% save the above linear discharge
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
%% calculate total sediment flux (interpolated and measured)
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

