% Compare mid-channel grain size data for BR

clear all,close all,clc
load('BR_TransectData_bed.mat')

% define the right and left bank coordinates of your idealized transect in UTM
[right.x,right.y]=deg2utm(16.099239, 95.320838);
[left.x,left.y]=deg2utm(16.100713, 95.330033);
deg2rotate=90+atand(abs(left.x-right.x)/abs(left.y-right.y));

for qq = 1:length(brt)
    [xx, yy,~] = deg2utm(brt(qq).lat,brt(qq).long);

    prj=proj([left.x - xx, left.y - yy],...
        [left.x- right.x, left.y - right.y]);
    
    brt(qq).dist = sqrt(prj(1).^2 + prj(2).^2);
    if brt(qq).dist<250
        brt(qq).stn=1;
    elseif brt(qq).dist>250 && brt(qq).dist<650
        brt(qq).stn=2;
    elseif brt(qq).dist>650
        brt(qq).stn=3;
    end
end
hf = brt(1:24);
lf = brt(25:end);
stn_hf=vertcat(hf.stn);
stn_lf=vertcat(lf.stn);

[stn_hf,srt] = sort(stn_hf,'ascend');
hf=hf(srt);
[stn_lf,srt] = sort(stn_lf,'ascend');
lf=lf(srt);

dist_hf=vertcat(hf.dist);
dist_lf=vertcat(lf.dist);

css_hf=vertcat(hf.frac);
css_lf=vertcat(lf.frac);
figure;
subplot(211),bar(dist_hf,css_hf,5,'stacked','EdgeColor',[1 1 1])
xlim([50 215])
subplot(212),bar(dist_lf,css_lf,5,'stacked','EdgeColor',[1 1 1])
xlim([50 215])
legend('clay','silt','sand')

%% Pattern on muddy flank:
flk_hf=mean(vertcat(hf(stn_hf==1 | stn_hf==3).frac))
flk_lf=mean(vertcat(lf(stn_lf==1 | stn_lf==3).frac))

%%
% d50_hf=vertcat(hf.median);
% d50_lf=vertcat(lf.median);
% figure;
% plot(stn_hf,d50_hf,'ko'),hold on
% plot(stn_lf,d50_lf,'ro')
% 
% figure;
% plot(dist_hf,d50_hf,'ko'),hold on
% plot(dist_lf,d50_lf,'ro')
