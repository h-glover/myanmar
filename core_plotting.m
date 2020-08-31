%% check the results by plotting results:
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

figure;
for jj=1:length(cores)
    subplot(1,6,jj)
    semilogx(act{jj},depth{jj},'o')
    ax=gca; ax.YDir='reverse';

    xlim([0.1 5])
    title(cores{jj})
end

subplot(163),xlabel('Excess activity (dpm/g), supp=0.7')
subplot(161),ylabel('Depth (cm)')