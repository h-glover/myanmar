clear all,close all,clc

load('C:\GLOVER\output\myanmar\grainsize\BR_GrainSize_Cores2.mat')
gs = cores.br463;
cd C:\GLOVER\output\myanmar\pb210
load('BR463_pb210ex7.mat')

figure;
depth = vertcat(pb210.midinterval);
act = vertcat(pb210.actatcoll_saltmudcor);
act(act<0)=NaN;
act_err = vertcat(pb210.tot_err)./2;
depth_err = vertcat(pb210.intervalsize)./2;

subplot(221)
errorbar(act,depth,depth_err,depth_err,act_err,act_err,'.')
axis ij
ax=gca; ax.XScale='log';
xlim([0.05 5]),ylim([0 100])

d_top = [5]; 
d_bot = [85];

% Code to define the data range over which to fit the line 
xx=act(depth>d_top & depth<d_bot & ~isnan(act));
yy=depth(depth>d_top & depth<d_bot &~isnan(act));
% Linear Regression
b=polyfit(log10(xx),yy,1);
x_fit=[min(xx),max(xx)];
% if x_fit(1)<0
%     x_fit(1)=0;
% end
y_fit = polyval(b,log10(x_fit));

hold on,
% plot(xx,yy,'ro')
semilogx(x_fit,y_fit,'k-')
title('Meinmahla Island')

% Text
[~,maxy]=max(y_fit);
[~,miny]=min(y_fit);
s =(0.0311*range(y_fit))/(log(x_fit(2)/x_fit(1)));
s = round(s,2)

y_fit = polyval(b,log10(xx));
r = corrcoef(yy,y_fit);
R2 = r(2,1)^2

% plot grain size
frac = vertcat(gs.frac)/100;
d50 = mean(vertcat(gs.median));

subplot(222)
barh(depth,frac,'stacked')
axis ij
ylim([0 100]), xlim([0 1])

%% load the second core and do all the same:
clearvars -except cores
gs = cores.br524;
gs([1,3,5])=[];
cd C:\GLOVER\output\myanmar\pb210
load('BR524_pb210ex7.mat')

depth = vertcat(pb210.midinterval);
act = vertcat(pb210.actatcoll_saltmudcor);
act(act<0)=NaN;
act_err = vertcat(pb210.tot_err)./2;
depth_err = vertcat(pb210.intervalsize)./2;

subplot(223)
errorbar(act,depth,depth_err,depth_err,act_err,act_err,'.')
axis ij
ax=gca; ax.XScale='log';
xlim([0.05 5]),ylim([0 100])

d_top = [0]; 
d_bot = [85];

% Code to define the data range over which to fit the line 
xx=act(depth>d_top & depth<d_bot & ~isnan(act));
yy=depth(depth>d_top & depth<d_bot &~isnan(act));
% Linear Regression
b=polyfit(log10(xx),yy,1);
x_fit=[min(xx),max(xx)];
% if x_fit(1)<0
%     x_fit(1)=0;
% end
y_fit = polyval(b,log10(x_fit));

hold on,
% plot(xx,yy,'ro')
% semilogx(x_fit,y_fit,'k-')
title('Agricultural Field')
% Text
[~,maxy]=max(y_fit);
[~,miny]=min(y_fit);
s =(0.0311*range(y_fit))/(log(x_fit(2)/x_fit(1)));
s = round(s,2)

y_fit = polyval(b,log10(xx));
r = corrcoef(yy,y_fit);
R2 = r(2,1)^2

% plot grain size
frac = vertcat(gs.frac)/100;
d50 = mean(vertcat(gs.median));

subplot(224)
barh(depth,frac,'stacked')
axis ij
ylim([0 100]), xlim([0 1])

legend({'clay','silt','sand'})

