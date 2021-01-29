clear all,close all,clc
load('MMsm_Sept17_aqd')
smch=aqd;
load('MMLG_Sept17_aqd')
lgch=aqd; clear aqd

% smooth water level using convolution
nn=3;
cc=[1/nn;1/nn;1/nn];
smch.pres(nn:end-nn)=conv(smch.pres(nn:end-nn),cc,'same');
lgch.pres(nn:end-nn)=conv(lgch.pres(nn:end-nn),cc,'same');

% calculate slope of pressure in m/hr
lgch.slope=lgch.pres(2:end)-lgch.pres(1:end-1);
lgch.slope=[lgch.slope(1);lgch.slope]*60;
lgch.slope(abs(lgch.slope)>1)=NaN;
lgch.slope(nn:end-nn)=conv(lgch.slope(nn:end-nn),cc,'same');
smch.slope=smch.pres(2:end)-smch.pres(1:end-1);
smch.slope=[smch.slope(1);smch.slope]*60;
smch.slope(abs(smch.slope)>1.5)=NaN;
smch.slope(nn:end-nn)=conv(smch.slope(nn:end-nn),cc,'same');

% figure;
% subplot(211)
% yyaxis left
% plot(lgch.time,lgch.pres,'k:'),hold on
% yyaxis right
% plot(lgch.time,lgch.slope,'k-')
% datetick('x','dd','keeplimits')
% 
% subplot(212)
% yyaxis left
% plot(smch.time,smch.pres,'k:'),hold on
% yyaxis right
% plot(smch.time,smch.slope,'k-')
% datetick('x','dd','keeplimits')

t_cycles=[216:(12.4*60):(6*12.4*60)];
for jj=1:length(t_cycles)-1
    V=[];
    V=lgch.slope(t_cycles(jj):t_cycles(jj+1));
    V=V(~isnan(V));
    mu3(jj)=(1/length(V))*sum((V.^3));
end
V=lgch.slope(~isnan(lgch.slope));
mu2=(1/length(V))*sum((V.^2));
lgch.gamma0=mu3./(mu2^(2/3));
lgch.gammatime=lgch.time(t_cycles(1:end-1));

t_cycles=[1:(12.4*60):(6*12.4*60)];
for jj=1:length(t_cycles)-1
    V=[];
    V=smch.slope(t_cycles(jj):t_cycles(jj+1));
    V=V(~isnan(V));
    mu3(jj)=(1/length(V))*sum((V.^3));
end
V=smch.slope(~isnan(smch.slope));
mu2=(1/length(V))*sum((V.^2));
smch.gamma0=mu3./(mu2^(2/3));
smch.gammatime=smch.time(t_cycles(1:end-1));


figure;
%subplot(211)
yyaxis left
plot(lgch.time(3:end-3),lgch.pres(3:end-3),'b:'),hold on
plot(smch.time(3:end-3),smch.pres(3:end-3),'k:')
yyaxis right
plot(lgch.gammatime,lgch.gamma0,'bo'),hold on
plot(smch.gammatime,smch.gamma0,'ko')


clear mu3
t_cycles=[1:(24.8*60):(4*24.8*60)];t_cycles(end)=length(lgch.time);
for jj=1:length(t_cycles)-1
    V=[];
    V=lgch.slope(t_cycles(jj):t_cycles(jj+1));
    V=V(~isnan(V));
    mu3(jj)=(1/length(V))*sum((V.^3));
end
lgch.gamma0_24=mu3./(mu2^(2/3));
lgch.gammatime_24=lgch.time(t_cycles(1:end-1));

t_cycles=[1:(24.8*60):(4*24.8*60)];t_cycles(end)=length(smch.time);
for jj=1:length(t_cycles)-1
    V=[];
    V=smch.slope(t_cycles(jj):t_cycles(jj+1));
    V=V(~isnan(V));
    mu3(jj)=(1/length(V))*sum((V.^3));
end
smch.gamma0_24=mu3./(mu2^(2/3));
smch.gammatime_24=smch.time(t_cycles(1:end-1));


plot(lgch.gammatime_24,lgch.gamma0_24,'bd')
plot(smch.gammatime_24,smch.gamma0_24,'kd')
ylim([0 0.5])
datetick('x','dd','keeplimits')
tt=round(min([lgch.time(1),smch.time(1)]));
%%
clearvars -except tt%,close all,clc
load('MMsm_Mar18_aqd')
smch=aqd;
load('MMLG_Mar18_aqd')
lgch=aqd; clear aqd

smch.time=smch.time-(round(smch.time(1))-tt);
lgch.time=lgch.time-(round(lgch.time(1))-tt);


% smooth water level using convolution
nn=3;
cc=[1/nn;1/nn;1/nn];
smch.pres(nn:end-nn)=conv(smch.pres(nn:end-nn),cc,'same');
lgch.pres(nn:end-nn)=conv(lgch.pres(nn:end-nn),cc,'same');

% calculate slope of pressure in m/hr
lgch.slope=lgch.pres(2:end)-lgch.pres(1:end-1);
lgch.slope=[lgch.slope(1);lgch.slope]*12;
lgch.slope(abs(lgch.slope)>1)=NaN;
lgch.slope(nn:end-nn)=conv(lgch.slope(nn:end-nn),cc,'same');
smch.slope=smch.pres(2:end)-smch.pres(1:end-1);
smch.slope=[smch.slope(1);smch.slope]*12;
smch.slope(abs(smch.slope)>1.5)=NaN;
smch.slope(nn:end-nn)=conv(smch.slope(nn:end-nn),cc,'same');

% figure;
% subplot(211)
% yyaxis left
% plot(lgch.time,lgch.pres,'k:'),hold on
% yyaxis right
% plot(lgch.time,lgch.slope,'k-')
% datetick('x','dd','keeplimits')
% 
% subplot(212)
% yyaxis left
% plot(smch.time,smch.pres,'k:'),hold on
% yyaxis right
% plot(smch.time,smch.slope,'k-')
% datetick('x','dd','keeplimits')
%
t_cycles=[1:round(12.4*12):round(4*12.4*12)];
for jj=1:length(t_cycles)-1
    V=[];
    V=lgch.slope(t_cycles(jj):t_cycles(jj+1));
    V=V(~isnan(V));
    mu3(jj)=(1/length(V))*sum((V.^3));
end
V=lgch.slope(~isnan(lgch.slope));
mu2=(1/length(V))*sum((V.^2));
lgch.gamma0=mu3./(mu2^(2/3));
lgch.gammatime=lgch.time(t_cycles(1:end-1));

t_cycles=[1:round(12.4*12):round(8*12.4*12)];
for jj=1:length(t_cycles)-1
    V=[];
    V=smch.slope(t_cycles(jj):t_cycles(jj+1));
    V=V(~isnan(V));
    mu3(jj)=(1/length(V))*sum((V.^3));
end
V=smch.slope(~isnan(smch.slope));
mu2=(1/length(V))*sum((V.^2));
smch.gamma0=mu3./(mu2^(2/3));
smch.gammatime=smch.time(t_cycles(1:end-1));


%subplot(212)
yyaxis left
plot(lgch.time(3:end-3),lgch.pres(3:end-3),'b--'),hold on
plot(smch.time(3:end-3),smch.pres(3:end-3),'k--')
yyaxis right
plot(lgch.gammatime,lgch.gamma0,'bo','MarkerFaceColor','b'),hold on
plot(smch.gammatime,smch.gamma0,'ko','MarkerFaceColor','k')

%
clear mu3
t_cycles=[1:round(24.8*12):round(2*24.8*12)];
for jj=1:length(t_cycles)-1
    V=[];
    V=lgch.slope(t_cycles(jj):t_cycles(jj+1));
    V=V(~isnan(V));
    mu3(jj)=(1/length(V))*sum((V.^3));
end
lgch.gamma0_24=mu3./(mu2^(2/3));
lgch.gammatime_24=lgch.time(t_cycles(1:end-1));

t_cycles=[1:round(24.8*12):round(5*24.8*12)];
t_cycles(end)=length(lgch.time);
for jj=1:length(t_cycles)-1
    V=[];
    V=smch.slope(t_cycles(jj):t_cycles(jj+1));
    V=V(~isnan(V));
    mu3(jj)=(1/length(V))*sum((V.^3));
end
smch.gamma0_24=mu3./(mu2^(2/3));
smch.gammatime_24=smch.time(t_cycles(1:end-1));


plot(lgch.gammatime_24,lgch.gamma0_24,'bd','MarkerFaceColor','b')
plot(smch.gammatime_24,smch.gamma0_24,'kd','MarkerFaceColor','k')
ylim([0 0.5])
datetick('x','dd','keeplimits')

