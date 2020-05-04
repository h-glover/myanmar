% Make plots of slope vs turbidity
clear all,close all,clc
F=dir('AQD*.mat');
cinv=repmat(1/3,3,1);
cc=1;
figure;
for jj=1:length(F)
load(F(jj).name)

% remove spikes
aqd.SSC1=medfilt1(aqd.SSC1,3);
aqd.SSC2=medfilt1(aqd.SSC2,3);


% calculate water slope in m/hr from depth
deltaT=round((1/24)./(aqd.time(2)-aqd.time(1)));
aqd.slope=(aqd.depth(2:end)-aqd.depth(1:end-1)).*deltaT;
aqd.slope(end+1)=aqd.slope(end);
aqd.slope(abs(aqd.slope)>3)=NaN;
aqd.slope(2:end-2)=conv(aqd.slope(2:end-2),cinv,'same');

subplot(3,3,cc)
scatter(aqd.slope,aqd.SSC1./1000,'k.')
hold on,scatter(aqd.slope,aqd.SSC2./1000,'r.')
axis tight
xlim([-1.5 1.5])
title(F(jj).name([5:9,11:12]))

if jj==3
    cc=cc+2;
else
    cc=cc+1;
end

end
subplot(3,3,1)
ylabel('SSC g/L')
xlabel('slope m/hr')

%% Make High and Low FLow small channel plots only, omitting data >2g/L
clear all,close all,clc
F=dir('AQD*.mat');

% remove LG and organize:
F=F([3,1,5,6]);
clr=[0.7 0.1 0; 0.7 0.1 0; 0.1 0.1 0.4; 0.1 0.1 0.4];% change color by season

% load and plot the low flow vals
for jj=1:2
load(F(jj).name)

% remove spikes
dd(jj).SSC1=medfilt1(aqd.SSC1,5);
dd(jj).SSC2=medfilt1(aqd.SSC2,5);
% remove high vals
dd(jj).SSC1(dd(jj).SSC1>2000 | dd(jj).SSC1<0)=NaN;
dd(jj).SSC2(dd(jj).SSC2>2000 | dd(jj).SSC2<0)=NaN;

% calculate water slope in m/hr from depth
deltaT=round((1/24)./(aqd.time(2)-aqd.time(1)));
% aqd.depth(2:end-2)=conv(aqd.depth(2:end-2),cinv,'same');
dd(jj).slope=(aqd.depth(2:end)-aqd.depth(1:end-1)).*deltaT;
dd(jj).slope(end+1)=dd(jj).slope(end);
dd(jj).slope=medfilt1(dd(jj).slope,5);
dd(jj).slope(abs(dd(jj).slope)>3)=NaN;


subplot(1,2,jj)
scatter(dd(jj).slope,dd(jj).SSC1./1000,...
    'Marker','d','MarkerEdgeColor',clr(jj,:))
hold on
scatter(dd(jj).slope,dd(jj).SSC2./1000,...
    'Marker','d','MarkerEdgeColor',clr(jj,:))
end


cinv=repmat(1/11,1,11);
for jj=3:4
    load(F(jj).name)
% remove spikes
dd(jj).SSC1=medfilt1(aqd.SSC1,5);
dd(jj).SSC2=medfilt1(aqd.SSC2,5);
% remove high vals
dd(jj).SSC1(dd(jj).SSC1>2000 | dd(jj).SSC1<0)=NaN;
dd(jj).SSC2(dd(jj).SSC2>2000 | dd(jj).SSC2<0)=NaN;

% calculate water slope in m/hr from depth
deltaT=round((1/24)./(aqd.time(2)-aqd.time(1)));
aqd.depth(2:end-2)=conv(aqd.depth(2:end-2),cinv,'same');
dd(jj).slope=(aqd.depth(2:end)-aqd.depth(1:end-1)).*deltaT;
dd(jj).slope(end+1)=dd(jj).slope(end);
dd(jj).slope=medfilt1(dd(jj).slope,5);
dd(jj).slope(abs(dd(jj).slope)>3)=NaN;
    
subplot(1,2,jj-2)
scatter(dd(jj).slope,dd(jj).SSC1./1000,...
    'Marker','.','MarkerEdgeColor',clr(jj,:))
hold on
scatter(dd(jj).slope,dd(jj).SSC2./1000,...
    'Marker','.','MarkerEdgeColor',clr(jj,:))
axis([-1.5 1.5 0 2])
end

