% silt into monkey pagoda interior:
% MP to dock = 346, 521, 347, 520, 522
clear all,close all,clc

cd C:\GLOVER\output\myanmar\grainsize
load('BR_GrainSize_Cores.mat')
figure;

% subplot(411)
% L = length(cores.br521);
% C = cmocean('amp',L);
% for jj=1:L
%     semilogx(cores.br521(jj).data(:,1),...
%         cores.br521(jj).data(:,2),'Color',C(jj,:)),hold on
% end
crs = {'br521','br347','br520','br522'};

for jj=1:4
    subplot(4,1,jj)
    frac = vertcat(cores.(crs{jj}).frac);
    lbl = vertcat(cores.(crs{jj}).name);
    lbl = lbl(:,7:11);
    dpth = (str2num(lbl(:,4:5)) - str2num(lbl(:,1:2)))/2 + str2num(lbl(:,1:2));
    bar(dpth,frac)
    ax=gca;ax.XTick = dpth; ax.XTickLabels = lbl;
    xlim([-0.5 70])
    ylim([0 65]),ylabel('Perc')

end
legend({'clay','silt','sand'}),title('Near Dock')
subplot(411),title('Pagoda')

