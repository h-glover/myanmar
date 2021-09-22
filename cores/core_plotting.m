%% check the results by plotting results:
clear all,close all,clc

cd C:\GLOVER\output\myanmar\pb210
F=dir('*ex.mat');
% F = F([1:5]);

for jj=1:length(F)
    load(F(jj).name)
    depth{jj} = vertcat(pb210.midinterval);
    act{jj} = vertcat(pb210.actatcoll_saltmudcor);
    act{jj}(act{jj}(:)<0.005)=NaN;
    cores{jj} = F(jj).name(1:5);
end

% figure;
for jj=1:length(cores)
    figure;
%     subplot(2,8,jj)
    semilogx(act{jj},depth{jj},'o')
    ax=gca; ax.YDir='reverse';
    axis tight, %ylim([0 100]),xlim([0.1 3])%,ylim([0 100])
    title(cores{jj})
end

% subplot(261)
% xlabel('Total activity (dpm/g)')
% ylabel('Depth (cm)')

%% plot the raw data without excess activity (supp=0)
clear all,close all,clc

cd C:\GLOVER\output\myanmar\pb210
F=dir('*ex.mat');

for jj=1:length(F)
    load(F(jj).name)
    depth{jj} = vertcat(pb210.midinterval);
    act{jj} = vertcat(pb210.actatcoll_saltmudcor);
    act{jj}(act{jj}(:)<0.001)=NaN;
    cores{jj} = F(jj).name(1:5);
end


for jj=1:length(cores)
    figure;%subplot(1,6,jj)
    semilogx(act{jj},depth{jj},'o')
    ax=gca; ax.YDir='reverse';

    xlim([0.1 5])
    title(cores{jj})
end

% subplot(163),xlabel('Total Activity (dpm/g)')
% subplot(161),ylabel('Depth (cm)')


%% Error bar plot of 5 best:
clear all,close all,clc

load('C:\GLOVER\output\myanmar\grainsize\BR_GrainSize_Cores2.mat')
cores.br346 = cores.br346([1,3:9]);
cores.br427 = cores.br427([1,2,4,6:8]);
cores.br520 = cores.br520([2:8]);
cores.br521 = cores.br521([1,2,4:10]);
cores.br522 = cores.br522([2:9]);
corename = {'br346','br427','br520','br521','br522'};

cd C:\GLOVER\output\myanmar\pb210
F=dir('*ex.mat');
F = F([1,4,10:12]);

figure;
for jj=1:5
load(F(jj).name)

depth = vertcat(pb210.midinterval);
act = vertcat(pb210.actatcoll_saltmudcor);
act(act<0)=NaN;
act_err = vertcat(pb210.tot_err)./2;
depth_err = vertcat(pb210.intervalsize)./2;


subplot(2,5,jj)
errorbar(act,depth,depth_err,depth_err,act_err,act_err,'.')
axis ij
ax=gca; ax.XScale='log';
xlim([0.1 5]),ylim([0 100])

d_top = [0,7,5,0,0]; 
d_bot = [30,50,30,20,40];

% Code to define the data range over which to fit the line 
xx=act(depth>d_top(jj) & depth<d_bot(jj) & ~isnan(act));
yy=depth(depth>d_top(jj) & depth<d_bot(jj) &~isnan(act));
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

% Text
[~,maxy]=max(y_fit);
[~,miny]=min(y_fit);
s =(0.0311*range(y_fit))/(log(x_fit(2)/x_fit(1)));
s = round(s,2);
title([F(jj).name(1:5),', A=',num2str(s),' cm/y'])


% plot grain size
frac = cores.(corename{jj});
frac = vertcat(frac.frac)/100;

subplot(2,5,jj+5)
barh(depth,frac,'stacked')
axis ij
ylim([0 100]), xlim([0 1])

end
colormap('gray')
legend({'clay','silt','sand'})


%% plots that dont have a log-linear section
clear all,close all,clc

load('C:\GLOVER\output\myanmar\grainsize\BR_GrainSize_Cores2.mat')
corename = {'br360','br462','br463'};

cd C:\GLOVER\output\myanmar\pb210
F=dir('*ex.mat');
F = F([3,8,9]);

figure;
for jj=1:length(F)
load(F(jj).name)

depth = vertcat(pb210.midinterval);
act = vertcat(pb210.actatcoll_saltmudcor);
act(act<0)=NaN;
act_err = vertcat(pb210.tot_err)./2;
depth_err = vertcat(pb210.intervalsize)./2;
 
subplot(2,3,jj)
errorbar(act,depth,depth_err,depth_err,act_err,act_err,'.')
axis ij
ax=gca; ax.XScale='log';
xlim([0.1 5]),ylim([0 100])

title(F(jj).name(1:5))

% plot grain size
frac = cores.(corename{jj});
frac = vertcat(frac.frac)/100;

subplot(2,3,jj+3)
barh(depth,frac,'stacked')
axis ij
ylim([0 100]), xlim([0 1])
end

%% plot the accum rate without grain size just for testing:

clear all,close all,clc

cd C:\GLOVER\output\myanmar\pb210
% 428, 432

F=dir('*ex.mat');
F = F([1,4,5,8,12:14]);
d_top = [0,7,0,0,5,0,0]; 
d_bot = [30,50,60,60,30,20,40];
% d_top = [0,7,5,0,0]; 
% d_bot = [30,50,30,20,40];


figure;
for jj=1:length(F)
load(F(jj).name)

depth = vertcat(pb210.midinterval);
act = vertcat(pb210.actatcoll_saltmudcor);
act(act<0)=NaN;
act_err = vertcat(pb210.tot_err)./2;
depth_err = vertcat(pb210.intervalsize)./2;

if jj==3
    act(3)=NaN;
end

subplot(2,4,jj)
errorbar(act,depth,depth_err,depth_err,act_err,act_err,'.')
axis ij
ax=gca; ax.XScale='log';
xlim([0.1 5]),ylim([0 100])


% Code to define the data range over which to fit the line 
xx=act(depth>d_top(jj) & depth<d_bot(jj) & ~isnan(act));
yy=depth(depth>d_top(jj) & depth<d_bot(jj) &~isnan(act));
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

% Text
[~,maxy]=max(y_fit);
[~,miny]=min(y_fit);
s =(0.0311*range(y_fit))/(log(x_fit(2)/x_fit(1)));
s = round(s,2);
title([F(jj).name(1:5),', A=',num2str(s),' cm/y'])


y_fit = polyval(b,log10(xx));
r = corrcoef(yy,y_fit);
R2(jj) = r(2,1)^2;

end
