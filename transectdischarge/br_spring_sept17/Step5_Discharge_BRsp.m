clear all
load BR_Sept17_spring.mat
for ii=1:length(adcp)
MeasuredDis(ii)=nansum(nansum(adcp(ii).alongComplete./100)).*0.5;

MeasuredTime(ii)=mean(time);
end

InterpTimes=days([min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.4)]);
InterpDis=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],[MeasuredDis,MeasuredDis(1)],InterpTimes,'linear');
InterpRes=mean(InterpDis);

save('BR_spring_discharge_interp','InterpDis','InterpTimes','InterpRes')
%%
clear all
load BR_Sept17_spring_SedFlux.mat
for ii=1:length(adcp)
MeasuredTime(ii)=mean(adcp(ii).time);
end
MeasuredDis=horzcat(adcp.MeasuredDis);
InterpTimes=days([min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.4)]);
InterpDis=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],[MeasuredDis,MeasuredDis(1)],InterpTimes,'linear');
InterpRes=mean(InterpDis);

MeasuredFlux=horzcat(adcp.sigmaMeasuredSedFlux);
InterpFlux=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],[MeasuredFlux,MeasuredFlux(1)],InterpTimes,'linear');
InterpResFlux=mean(InterpFlux);


save('BR_spring_discharge_interp','Interp*','Measured*')
