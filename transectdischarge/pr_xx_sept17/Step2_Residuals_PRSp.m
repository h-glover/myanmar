clear all
% load the mat file containng the adcp structure produced by running Dan Nowacki's adcp_ascii function 
load PR_Sept17_Spring.mat
% PR transect locations % fix the ordering of the transects:
% NPD transect: (1)
% PR_long_002,003 = adcp(1:2)
% 
% Main Transect:(2)
% PR_long2_002,003
% xPR_trans_001,002,005,006
% = adcp(3:7)
% 
% PR, above split Transect:(3)
% xPR_trans2_002,003
% =adcp(8:9)

c=1;
UpDists=[];
DownDists=[];
AllVels=[];
transect_length=2200; % in meters
allalong=NaN(max(vertcat(adcp.numcells)),transect_length+1,48);
allacross=NaN(max(vertcat(adcp.numcells)),transect_length+1,48);


% define the right and left bank coordinates of your idealized transect in UTM 
[right(1).x,right(1).y]=deg2utm(16.5487, 94.6834);% Near NPD
[left(1).x,left(1).y]=deg2utm(16.5444,94.6894);% Near NPD
[right(2).x,right(2).y]=deg2utm(16.48402, 94.67073);% Main Trans
[left(2).x,left(2).y]=deg2utm(16.48433, 94.69185);% Main Trans
[right(3).x,right(3).y]=deg2utm(16.50069, 94.66892);% PR, south of NPD
[left(3).x,left(3).y]=deg2utm(16.50143, 94.68117);% PR, south of NPD

for jj=1:length(right)
deg2rotate(jj)=90+atand(abs(left(jj).x-right(jj).x)/abs(left(jj).y-right(jj).y));
end
%
for n=1:length(adcp)
  
    % convert the coordinates of of your ADCP records to UTM
    [adcp(n).x, adcp(n).y,~] = deg2utm(adcp(n).lat(~isnan(adcp(n).lat)), adcp(n).lon(~isnan(adcp(n).lon)));


    % select correct transect for PR
    if n<=2
        jj=1;
    elseif n>2 && n<=7
        jj=2;
    elseif n>7
        jj=3;
    end
    % project the ADCP data onto the idealized transect
    for qq = 1:length(adcp(n).x)
        adcp(n).proj(qq,:) = proj([left(jj).x - adcp(n).x(qq),...
            left(jj).y - adcp(n).y(qq)], [left(jj).x- right(jj).x,...
            left(jj).y - right(jj).y]);
    end

    % define new distances along the idealized transect
    adcp(n).dist =  sqrt(adcp(n).proj(:,1).^2 + adcp(n).proj(:,2).^2);
    
    % perform the coordinate rotation into along and accross components
    [adcp(n).across, adcp(n).along] = rot_earth(adcp(n).east, adcp(n).north, deg2rotate(jj));
    
    % define transect width 
    adcp(n).interpdist = [0:1:transect_length]';
    
    % find only the unique data (eliminate distance repeats)
    [adcp(n).uniquedist,uniqueindices]=unique(adcp(n).dist);
    
    %find unique depths 
    
    adcp(n).meandepth=mean(adcp(n).depth);
    adcp(n).uniquedepths=adcp(n).meandepth(uniqueindices);
    
    % interpolate the 
    adcp(n).interpalong = interp1(adcp(n).uniquedist, adcp(n).along(:,uniqueindices)', adcp(n).interpdist)';
    adcp(n).interpacross = interp1(adcp(n).uniquedist, adcp(n).across(:,uniqueindices)', adcp(n).interpdist)';
    adcp(n).interpdepths = interp1(adcp(n).uniquedist, adcp(n).meandepth(:,uniqueindices)', adcp(n).interpdist)';
    % compute x-sectional area
    
     
    c=c+1;
end
save PR_Sept17_Spring.mat adcp

%%

for i=1:length(adcp)
figure(100+i)
clf 
set(gcf,'color','w')
% p=pcolor(adcp(n).interpdist,adcp(n).z(:,1),meanalong);
% if i==15
%     p=surf(adcp(i).interpdist,adcp(2).z(:,1),(smooth2a(allalong(:,:,i).*-1,4,10)));
%     shading flat
%     view(0,90)
%     c=colorbar;
%     caxis([-200,200])
%     set(gca,'ydir','reverse')
%     ylim([0,20])
%     xlim([0,850])
%     set(gca,'fontsize',textsize)
%     title('Yangon River September 2017','fontsize',textsize)
%     ylabel(c, 'Downstream Velocity (cm s^-^1)','fontsize',textsize);
%     xlabel('Cross-Channel Distance (m)','fontsize',textsize);
%     ylabel('Water Depth (m)','fontsize',textsize);
%     colormap(othercolor('BrBG8'))
%     Frames(i)=getframe(gcf);
% else
p=surf(adcp(i).interpdist,adcp(1).z(:,1),(smooth2a(allalong(:,:,i),4,4)));
% set(gca,'ydir','reverse')
shading flat
view(0,90)
c=colorbar;
caxis([-200,300])
set(gca,'ydir','reverse')
ylim([0,20])
xlim([0,850])
set(gca,'fontsize',textsize)
% title(FileList(i+2).name,'fontsize',textsize)
title('Yangon River | Spring Tide | High Discharge','fontsize',textsize)
ylabel(c, 'Downstream Velocity (cm s^-^1)','fontsize',textsize);
xlabel('Cross-Channel Distance (m)','fontsize',textsize);
ylabel('Water Depth (m)','fontsize',textsize);
colormap(othercolor('BrBG8'))
Frames(i)=getframe(gcf);
end

 %%
% 
%  % create the video writer with 1 fps
%  writerObj = VideoWriter('YangonRiverSept2017.avi');
%  writerObj.FrameRate = 1;
%  % set the seconds per image
%  secsPerImage = ones(length(Frames)).*1;
%  % open the video writer
%  open(writerObj);
%  % write the frames to the video
%  for u=2:length(Frames)
%      for v=1:secsPerImage(u) 
%          writeVideo(writerObj, Frames(u));
%      end
%  end
%  % close the writer object
%  close(writerObj);


%%
% figure(4)
% clf 
% set(gcf,'color','w')
% p=pcolor(adcp(n).interpdist,adcp(2).z(:,1),smooth2a(meanalong,4,10));
% shading flat
% hold on
% cont=contour(adcp(n).interpdist,adcp(2).z(:,1),smooth2a(meanalong,4,20),[0,0],'linecolor','k','linewidth',.5);
% c=colorbar
% set(gca,'ydir','reverse')
% ylim([0,25])
% caxis([-200,200])
% set(gca,'fontsize',textsize)
% title('Yangon River - High Flow','fontsize',textsize)
% ylabel(c, 'Residual Downstream Velocity (cm s^-^1)','fontsize',textsize);
% xlabel('Cross-Channel Distance (m)','fontsize',textsize);
% ylabel('Water Depth (m)','fontsize',textsize);
% colormap(othercolor('BrBG8'))
% set(gcf,'Position',[100,600,800,400])
% text(200,23,'\leftarrow West','fontsize',textsize)
% text(7750,23,'East \rightarrow','fontsize',textsize)
% export_fig YangonRiverResidual.pdf -nocrop
% 
% figure(45)
% clf 
% set(gcf,'color','w')
% p=pcolor(adcp(n).interpdist,adcp(n).z(:,1),smooth2a(meanacross,4,6));
% shading flat
% hold on
% cont=contour(adcp(n).interpdist,adcp(n).z(:,1),smooth2a(meanacross,4,6),[0,0],'linecolor','k','linewidth',.5)
% c=colorbar
% set(gca,'ydir','reverse')
% ylim([0,35])
% clim([-10,10])
% set(gca,'fontsize',textsize)
% title('Xingu - Low Flow','fontsize',textsize)
% ylabel(c, 'Residual Across-Ría (east) Velocity (cm s^-^1)','fontsize',textsize);
% xlabel('Cross-Channel Distance (m)','fontsize',textsize);
% ylabel('Water Depth (m)','fontsize',textsize);
% colormap(othercolor('BrBG8'))
% set(gcf,'Position',[100,600,800,400])
% export_fig XinguResidual2011_across.pdf -nocrop

% figure(44)
% clf 
% set(gcf,'color','w')
% scatter(DownVels,UpVels)
% set(gca,'fontsize',textsize)
% xlabel('Downstream Mean Velocity','fontsize',textsize);
% ylabel('Upstream Mean Velocity','fontsize',textsize);
% 
% figure(46)
% clf 
% set(gcf,'color','w')
% plot(MeanVels)
% set(gca,'fontsize',textsize)
% xlabel('Crossings','fontsize',textsize);
% ylabel('Mean Velocity','fontsize',textsize);
% %%
% figure(47)
% clf
% set(gcf,'color','w')
% hist(AllVels,1000)
% xlim([-70,70])
% 
% h=histfit(AllVels);
% hist(UpDists,100)
% hold all
% hist(DownDists,100)

%%
% clear all




