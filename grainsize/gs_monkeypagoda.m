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
    
    % remove OM signal - shows up as sand
    for kk=1:length(cores.(crs{jj}))
        xx = cores.(crs{jj})(kk).frac;
        if xx(3)>9
            rr = xx(3)-9;
            xx(1)=xx(1) + (xx(1)/(xx(1)+xx(2))*rr);
            xx(2)=xx(2) + (xx(2)/(xx(1)+xx(2))*rr);
            xx(3) = 100 - (xx(1)+xx(2));
            cores.(crs{jj})(kk).frac = xx;
         end
    end
    frac = vertcat(cores.(crs{jj}).frac);
    
    % find the depth from the name:
    lbl = vertcat(cores.(crs{jj}).name);
    lbl = lbl(:,7:11);
    dpth = (str2num(lbl(:,4:5)) - str2num(lbl(:,1:2)))/2 + str2num(lbl(:,1:2));

    % plot
    bar(dpth,frac,'stacked')
    ax=gca;ax.XTick = dpth; ax.XTickLabels = lbl;
    xlim([-0.5 70])
    ylim([0 100]),ylabel('Perc')

end
legend({'clay','silt','sand'}),title('Near Dock')
subplot(411),title('Pagoda')

%% plot avgs in each depth interval:

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
% dpth_bin = [0,6,16,32,100];
% depth_bin_plot=[3,12,24,50];
% lbl_bin = {'0-6','6-16','16-32','>32'};
dpth_bin = [0,30,100];
lbl_bin = {'<32','>32'};
for jj=1:4
    subplot(4,1,jj)
    
    % remove OM signal - shows up as sand
    for kk=1:length(cores.(crs{jj}))
        xx = cores.(crs{jj})(kk).frac;
        if xx(3)>9
            rr = xx(3)-9;
            xx(1)=xx(1) + (xx(1)/(xx(1)+xx(2))*rr);
            xx(2)=xx(2) + (xx(2)/(xx(1)+xx(2))*rr);
            xx(3) = 100 - (xx(1)+xx(2));
            cores.(crs{jj})(kk).frac = xx;
         end
    end
    frac = vertcat(cores.(crs{jj}).frac);
    
    % find the depth from the name:
    lbl = vertcat(cores.(crs{jj}).name);
    lbl = lbl(:,7:11);
    dpth = (str2num(lbl(:,4:5)) - str2num(lbl(:,1:2)))/2 + str2num(lbl(:,1:2));
    
    % bin avg by depth in ~10-20 year intervals
    for kk=2:length(dpth_bin)
        idx = find(dpth>=dpth_bin(kk-1) & dpth<dpth_bin(kk));
        frac_bin(kk-1,:) = mean(frac(idx,:),1);
    end
    
    % plot
%     bar(depth_bin_plot,frac_bin,'stacked')
%     ax=gca;ax.XTick = depth_bin_plot; ax.XTickLabels = lbl_bin;
    bar(frac_bin,'stacked')
    ax=gca;ax.XTickLabels = lbl_bin;
    ylim([0 100]),ylabel('Perc')

end
legend({'clay','silt','sand'}),title('Near Dock')
subplot(411),title('Pagoda')


%% look at 326 core from near well:


clear all,close all,clc

cd C:\GLOVER\output\myanmar\grainsize
load('BR_GrainSize_Cores.mat')

cores.br335(3)=[];
figure;
frac = vertcat(cores.br335.frac);
lbl = vertcat(cores.br335.name);
lbl = lbl(:,7:11);
dpth = (str2num(lbl(:,4:5)) - str2num(lbl(:,1:2)))/2 + str2num(lbl(:,1:2));
% plot
    bar(dpth,frac,'stacked')
    ax=gca;ax.XTick = dpth; ax.XTickLabels = lbl;
    xlim([-0.5 70])
    ylim([0 100]),ylabel('Perc')

