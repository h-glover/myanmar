% calculae the shear stress over the tidal cycle using the law of the wall
% then the QSL
clear all,close all,clc

% first calculate the velocity profiles at the 3 locations (200, 488
% (thalweg), 800m)
load('BR_Sept17_Spring.mat')
adcp(4)=[];
% create a height-above-bed vector for the new bin coordinates (cm)
hab_mat=100:10:1800;
hab_log=log(hab_mat);
zdepth=round(adcp(1).z(:,1).*100);

for jj=1:length(adcp)
    % calculate the cross-section oriented speed(cm/s):
    adcp(jj).interpSpeed=sqrt(adcp(jj).interpalong.^2 + adcp(jj).interpacross.^2);
    
    % fix interp depth and convert to cm
    adcp(jj).interpdepths(...
        adcp(jj).interpdepths<0 | adcp(jj).interpdepths>26)=NaN;
    adcp(jj).interpdepths=round(adcp(jj).interpdepths.*100);

    % calculate the actual height above bed of each profile's bins (cm)
    adcp(jj).hab_actual=adcp(jj).interpdepths-zdepth;
    
    % for each column with depth>0, interpolate to the new HAB matrix so
    % that all of the profile bins line up
    cols=find(~isnan(adcp(jj).interpdepths));
    adcp(jj).spd_hab=NaN([length(hab_mat),length(adcp(jj).interpdepths)]);
    for nn=cols
        adcp(jj).spd_hab(:,nn)=interp1(...
            adcp(jj).hab_actual(:,nn),adcp(jj).interpSpeed(:,nn),hab_mat);
    end
    MeasuredTime(jj)=adcp(jj).time(1);
end

%% calculate ustar using a variety of depth ranges (using the whole profile
% through using half the profile) and save in the avgs structure
testavgs=[5:2:20];

for ii=1:length(testavgs)
xbins=5;%testavgs(ii);

for jj=1:length(adcp)
    % downsample the transect by averaging at 3 locations
    % within xbins of 200, within xbins of thalweg, within xbins of 800m
    spd_hab(:,1:3)=NaN(length(adcp(jj).spd_hab(:,1)),3);
    spd_hab(:,1)=nanmedian(adcp(jj).spd_hab(:,200-xbins:200+xbins),2);
    spd_hab(:,2)=nanmedian(adcp(jj).spd_hab(:,488-xbins:488+xbins),2);
    spd_hab(:,3)=nanmedian(adcp(jj).spd_hab(:,800-xbins:800+xbins),2);

    % find first nearbed bin with data
    [mm,~]=size(spd_hab);
    bin1=sum(isnan(spd_hab(1:20,:)))+1;
    binend=(mm-sum(isnan(spd_hab(20:end,:))))-testavgs(ii);
    binend(2)=binend(2)-5;
    
    % Calculate line fit to log(z) and ubar for the 3 stations (kk)
    
%     figure(jj)
    for kk=1:3
    xx=[];yy=[];vals=[];p=[];pfit=[];
    vals=[bin1(kk)+1,binend(kk)];
    xx=spd_hab(vals,kk);
    yy=hab_log(vals)';
%     subplot(1,3,kk)
%     plot(spd_hab(:,kk),hab_log,'k'),hold on
%     plot(xx,yy,'o')

    % calculate ubar to calculate Cd (cm/s)
    avgs(ii).ubar(jj,kk) = nanmean(spd_hab(bin1(kk):bin1(kk)+2,kk));
    avgs(ii).ubar_depth(jj,kk)=nanmean(hab_mat(bin1(kk):bin1(kk)+2));
    %plot(adcp(jj).ubar(kk),log(adcp(jj).ubar_depth(kk)),'*')
    
    % calculate the fit line
    pfit=polyfit(xx,yy,1);
    yfit=polyval(pfit,xx);
%     plot(xx,yfit,'r')

    % ustar is K/slope: cm/s
    avgs(ii).ustar(jj,kk)=0.41/pfit(1,1);
    % z0 is exp(intercept): cm
    avgs(ii).z0(jj,kk)=exp(pfit(1,2));
    % cd=ustar^2/ubar^2
    avgs(ii).Cd_lotw(jj,kk)=(avgs(ii).ustar(jj,kk).^2)/(avgs(ii).ubar(jj,kk).^2);

    % bed shear stress = rho*ustar^2=g/cm3*cm2/s2=g/cm/s2/10=kg/m/s2
    avgs(ii).tau_lotw(jj,kk)=(avgs(ii).ustar(jj,kk).^2)/10;

    
    if pfit(1,1)<0 || avgs(ii).tau_lotw(jj,kk)>12
    avgs(ii).ustar(jj,kk)=NaN;
    avgs(ii).z0(jj,kk)=NaN;
    avgs(ii).Cd_lotw(jj,kk)=NaN;
    avgs(ii).tau_lotw(jj,kk)=NaN;
    
    end
    end
    
    
end
    avgs(ii).time=MeasuredTime;
    %clearvars -except ii *avgs adcp hab*
end


% use Robin's method of calculating Cd from slope fitting between ustar-ubar

% ubar is the same for each calculation version
ubar=avgs(1).ubar;

for kk=1:length(avgs) % calc Cd for each calculation

%     figure;
for jj=1:3 % calc Cd for each station
    p=[];xx=[];yy=[];
    
    % fit a line through u*^2 and ubar^2
    xx=ubar(~isnan(avgs(kk).ustar(:,jj)),jj).^2;
    yy=avgs(kk).ustar(~isnan(avgs(kk).ustar(:,jj)),jj).^2;
   
    % slope of fit line is Cd
    p=polyfit(xx,yy,1);
    avgs(kk).Cd_qsl(jj)=p(1,1);
    
    % plot ustar-ubar and the fit line
%     pfit=polyval(p,xx);
%     subplot(1,3,jj)
%     scatter(xx,yy,'k'),hold on
%     plot(xx,pfit)

end
    %calculate tau using QSL
    avgs(kk).tau_qsl= avgs(kk).Cd_qsl(jj).*(ubar.^2)/10; 
end

save('BR_Sept17_Spring_ShearStress','avgs')

% take the avg of all the tau
tau_qsl=cat(3,avgs.tau_qsl);
tau_qsl_std=std(tau_qsl,0,3);
tau_qsl=nanmean(tau_qsl,3);
% plot the tau result

figure;
errorbar(tau_qsl,tau_qsl_std),legend({'200','thal','800'})
