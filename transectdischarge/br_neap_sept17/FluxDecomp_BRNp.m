% flux decomposition into Fluvial and Tidal components (simple)
clear all%,close all,clc
load('BR_Sept17_Neap_SedFlux.mat')

% 
for ii=1:length(adcp)
MeasuredTime(ii)=mean(adcp(ii).time);
end
InterpTimes=days([min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.4)]);

%% calculate total Area and tidally avgd area at the cross section
depths=vertcat(adcp.interpdepths);%get depth across section
depths(depths<0 | depths>20)=NaN;% remove bad vals and fill with nearest
depths=fillmissing(depths,'nearest',2);
A=sum(depths,2,'omitnan')'; %integrate: summ depths
% make a tidally uniform vector by interpolating to 1 min interval
A_interp=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],...
    [A,A(1)],InterpTimes,'linear');
A0=mean(A_interp);
%% calculate u0 from integrated velocity (Qr) at each time point from sigma coords
binheight=depths/1000; % convert depth to sigma coordinate bin height
for jj=1:length(adcp)
    %binheight(jj,:)=adcp(jj).sigmaDepth(3,:)-adcp(jj).sigmaDepth(2,:);
    u(jj)=nansum(nansum(adcp(jj).sigmaAlongComplete.*binheight(jj,:)));%sum(m/s*m2)=m3/s
end
u_interp=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],...
    [u,u(1)],InterpTimes,'linear');
u0=mean(u_interp)/A0;%m/s

%% calculate c0 from integrated SSC (Qr) at each time point from sigma coords
for jj=1:length(adcp)
    c(jj)=nansum(nansum(adcp(jj).ssc./1000.*binheight(jj,:)));%sum(kg/m3 * m2)=kg/m
end
c_interp=interp1([MeasuredTime,days(min(MeasuredTime)+hours(12.4))],...
    [c,c(1)],InterpTimes,'linear');
c0=mean(c_interp)/A0; %units of kg/m3


Fr=c0*u0*A0/1000 % kg/m3 * m/s * m2 = kg/s /1000 = t

%% plot the u0 c0 and A0

figure(1)
subplot(311),plot(u_interp),datetick('x','HH','keeplimits'),hold on
subplot(312),plot(c_interp),datetick('x','HH','keeplimits'),hold on
subplot(313),plot(A_interp),datetick('x','HH','keeplimits'),hold on

subplot(311),ylabel('int(u), m3/s')
subplot(312),ylabel('int(ssc), kg/m')
subplot(313),ylabel('Area, m^2')

%subplot(311),legend('HF Neap','HF Spring','LF Neap','LF Spring'),ylabel('int(u), m3/s')
