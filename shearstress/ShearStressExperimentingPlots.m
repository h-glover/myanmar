% experimental code and plotting

%% Compare drag coeffs 

clear all,close all,clc
F=dir('BR*.mat');
   figure;
for ff=1:4
   load(F(ff).name) 
   Cd(ff).lotw=vertcat(avgs.Cd_lotw);
   Cd(ff).qsl=vertcat(avgs.Cd_qsl);

   subplot(121),histogram(Cd(ff).lotw,[0:0.001:0.01]),hold on
   subplot(122),histogram(Cd(ff).qsl,[0:0.001:0.01]),hold on
end

%%
load('BR_Sept17_Spring.mat')
%adcp=adcp([1:8,35:end]); % ebb
%adcp=adcp([12:26]); % flood
%adcp([9:11,27:34])=[];% notrans
%adcp=adcp([9:11]); %ebbflood
%adcp=adcp([27:34]); %floodebb
%%
close all,clear all
load('shearstresscals.mat')
ubar=horzcat(avgs.ubar);
ustar=horzcat(avgs.ustar);

    %ubar=avgs(kk).ubar;
    %ustar=avgs(kk).ustar;
figure;
for jj=1:3
    p=[];xx=[];yy=[];
    idx=jj:3:length(ubar(1,:));
    xx=nanmean(ubar(:,idx),2);
    yy=nanmean(ustar(:,idx),2);
    xx=xx(~isnan(yy)).^2;
    yy=yy(~isnan(yy)).^2;
    p=polyfit(xx,yy,1);
    pfit=polyval(p,xx);
    subplot(1,3,jj)
    scatter(xx,yy,'k'),hold on
    plot(xx,pfit)
    Cd(jj)=p(1,1);
end

load('Ased.mat')

Ased(3).UmeanA350=100;
UstarA=vertcat(Ased.UStarA);
UmeanA350=vertcat(Ased.UmeanA350);

polyindex = find(~isnan(UstarA));
xx=(UmeanA350(polyindex)).^2;
yy=UstarA(polyindex).^2;
% Cd is slope from U*^2 against Umean^2
CdPoly350 = polyfit(xx,yy,1);
pfit=polyval(CdPoly350,xx);
figure;
scatter(xx,yy),hold on
plot(xx,pfit)

%%

clear all,close all,clc

load('BR_Sept17_Spring.mat')
load('DragCoeff')
load('BR_Sept17_Spring_FluxDecomp2.mat')

for jj=1:length(adcp)
    %add drag coefficient
    if jj<=8 || jj>=35
        adcp(jj).Cd=totalavg(:,2);
    elseif jj>=9 && jj<=11
        adcp(jj).Cd=totalavg(:,3);
    elseif jj>=12 && jj<=26
        adcp(jj).Cd=totalavg(:,4);
    elseif jj>=27 && jj<=34
        adcp(jj).Cd=totalavg(:,5);
    end
    MeasuredTime(jj)=adcp(jj).time(1);
end

Cd=horzcat(adcp.Cd);
spd=fluxdecomp.Qf./fluxdecomp.A;

C={'k','r','g'};
figure;
yyaxis left 
plot(fluxdecomp.time,spd)
yyaxis right
for jj=1:3
%subplot(3,1,jj)
plot(MeasuredTime,Cd(jj,:),'o','Color',C{jj}),hold on
end
ylabel('Cd'),ylim([0 0.01])
legend({'200','thal','800'})

%%
ustar=vertcat(avgs.ustar);
Cd=vertcat(avgs.Cd);
z0=vertcat(avgs.z0);

for ff=1:3
    figure(1)
    subplot(1,3,ff)
    histogram(Cd(:,ff),[0:0.001:0.01],...
        'DisplayStyle','stairs','normalization','probability')
    xlabel('cd'),hold on
    
    figure(2)
    subplot(1,3,ff)
    histogram(z0(:,ff),[0:1:10],...
        'DisplayStyle','stairs','normalization','probability')
    xlabel('z0'),hold on
    
    figure(3)
    subplot(1,3,ff)
    histogram(ustar(:,ff),'DisplayStyle','stairs','normalization','probability')
    xlabel('ustar'),hold on
end

%%
close all,clear all
F=dir('*cals*.mat');

for jj=[1,3]
    load(F(jj).name)

    ustar=vertcat(avgs.ustar);
Cd=vertcat(avgs.Cd);
z0=vertcat(avgs.z0);

for ff=1:3
    figure(1)
    subplot(1,3,ff)
    histogram(Cd(:,ff),[0:0.001:0.01],...
        'DisplayStyle','stairs','normalization','probability')
    xlabel('cd'),hold on
    
    figure(2)
    subplot(1,3,ff)
    histogram(z0(:,ff),[0:1:10],...
        'DisplayStyle','stairs','normalization','probability')
    xlabel('z0'),hold on
    
    figure(3)
    subplot(1,3,ff)
    histogram(ustar(:,ff),'DisplayStyle','stairs','normalization','probability')
    xlabel('ustar'),hold on
end

end
% 
% figure;
% for ff=1:3
% scatter(ustar(:,ff),z0(:,ff)),hold on
% end
% ylabel('z0'),xlabel('ustar')


%%
close all
ustar=horzcat(avgs.ustar);
Cd=horzcat(avgs.Cd_lotw); %Cd(Cd>0.01)=NaN;
z0=horzcat(avgs.z0);

for jj=1:3
    idx=jj:3:length(ustar(1,:));
    ustar_avg=nanmean(ustar(:,idx),2);
    Cd_avg=nanmean(Cd(:,idx),2);
    z0_avg=nanmean(z0(:,idx),2);
    
    figure(1)
    subplot(311),plot(ustar_avg),hold on
    subplot(312),plot(Cd_avg),hold on,%ylim([0 0.01])
    subplot(313),plot(z0_avg),hold on,%ylim([0 10])
    
    figure(1+jj)
    subplot(311),boxplot(ustar(:,idx)')
    subplot(312),boxplot(Cd(:,idx)'),%ylim([0 0.01])
    subplot(313),boxplot(z0(:,idx)'),%ylim([0 10])
end
    %figure(1),legend({'200','thal','800'})
%%
close all,clear all
F=dir('*cals*.mat');

for ff=[1,3]
    load(F(ff).name)
ustar=horzcat(avgs.ustar);
Cd=horzcat(avgs.Cd); %Cd(Cd>0.01)=NaN;
z0=horzcat(avgs.z0);

for jj=1:3
    idx=jj:3:length(ustar(1,:));
    ustar_avg=nanmean(ustar(:,idx),2);
    Cd_avg=nanmean(Cd(:,idx),2);
    z0_avg=nanmean(z0(:,idx),2);
    
    figure(ff)
    subplot(311),plot(ustar_avg),hold on
    subplot(312),plot(Cd_avg),hold on%,ylim([0 0.01])
    subplot(313),plot(z0_avg),hold on%,ylim([0 50])
    
    figure(ff*10+jj)
    subplot(311),boxplot(ustar(:,idx)')
    subplot(312),boxplot(Cd(:,idx)')%,ylim([0 0.01])
    subplot(313),boxplot(z0(:,idx)')%,ylim([0 50])
end
    figure(1),legend({'200','thal','800'})
end


% figure;
% plot(fluxdecomp.time,fluxdecomp.Qf)

%%
close all,clear all
F=dir('*cals*.mat');
cc=1;
for ff=[1:5]
    load(F(ff).name)
    ustar=horzcat(avgs.ustar);
    Cd=horzcat(avgs.Cd); Cd(ustar>7)=NaN;
    z0=horzcat(avgs.z0);

for jj=1:3
    idx=jj:3:length(ustar(1,:));
    Cd_avg=nanmedian(Cd(:,idx),2);
    totalavg(jj,ff)=nanmedian(Cd_avg);
    
    figure(10)
    subplot(5,1,cc),plot(Cd_avg),hold on,ylim([0 0.01])
    refline(0,totalavg(jj,ff))
    
    figure(ff)
    subplot(3,1,jj),boxplot(Cd(:,idx)'),ylim([0 0.01])

end
    %figure(10),subplot(311),legend({'200','thal','800'})
    cc=cc+1;
end
key={'all,1','ebb,1','ebbflood,1','flood,1','floodebb,1';...
    'all,2','ebb,2','ebbflood,2','flood,2','floodebb,2';...
    'all,3','ebb,3','ebbflood,3','flood,3','floodebb,3'};
% yyaxis right
% plot(MeasuredTime,1:length(MeasuredTime),'o')


%%
load('DragCoeff_Neap')

%adcp=adcp([11:26]); % ebb
%adcp=adcp([33:44]); % flood
%adcp([11:26,33:44])=[];% notrans
%adcp=adcp([26:32]); %ebbflood
%adcp=adcp([1:10,45:end]); %floodebb

% create a height-above-bed vector for the new bin coordinates (cm)
hab_mat=100:10:1800;
hab_log=log(hab_mat);
zdepth=round(adcp(1).z(:,1).*100);
xbins=5;

for jj=1:length(adcp)
    %add drag coefficient
    if jj>=11 && jj<=26
        adcp(jj).Cd=totalavg(:,2);
    elseif jj>=26 && jj<=32
        adcp(jj).Cd=totalavg(:,3);
    elseif jj>=33 && jj<=44
        adcp(jj).Cd=totalavg(:,4);
    elseif jj<=10 || jj>=45
        adcp(jj).Cd=totalavg(:,5);
    end
    
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

    % downsample the transect by averaging at 3 locations
    % within xbins of 200, within xbins of thalweg, within xbins of 800m
    spd_hab(:,1:3)=NaN(length(adcp(jj).spd_hab(:,1)),3);
    spd_hab(:,1)=nanmedian(adcp(jj).spd_hab(:,200-xbins:200+xbins),2);
    spd_hab(:,2)=nanmedian(adcp(jj).spd_hab(:,488-xbins:488+xbins),2);
    spd_hab(:,3)=nanmedian(adcp(jj).spd_hab(:,800-xbins:800+xbins),2);

    % find first nearbed bin with data
    [mm,~]=size(spd_hab);
    bin1=sum(isnan(spd_hab(1:20,:)))+1;
    binend=(mm-sum(isnan(spd_hab(20:end,:))))-15;
    
    % Calculate line fit to log(z) and ubar
    p=[];pfit=[];
%     figure(jj)
    for kk=1:3
    xx=[];yy=[];vals=[];
    vals=[bin1(kk)+1,binend(kk)];
    xx=spd_hab(vals,kk);
    yy=hab_log(vals)';
%     subplot(1,3,kk)
%     plot(spd_hab(:,kk),hab_log,'k'),hold on
%     plot(xx,yy,'o')

    % calculate ubar to calculate Cd (cm/s)
    adcp(jj).ubar(kk) = nanmean(spd_hab(bin1(kk):bin1(kk)+2,kk));
    adcp(jj).ubar_depth(kk)=nanmean(hab_mat(bin1(kk):bin1(kk)+2));
%     plot(adcp(jj).ubar(kk),log(adcp(jj).ubar_depth(kk)),'*')
    
    % calculate the fit line
    pfit(kk,1:2)=polyfit(xx,yy,1);
    yfit=polyval(pfit(kk,1:2),xx);
%     plot(xx,yfit,'r')

    % ustar is K/slope: cm/s
    adcp(jj).ustar(kk)=0.41/pfit(kk,1);

    % bed shear stress = rho*ustar^2=g/cm3*cm2/s2=g/cm/s2/10=kg/m/s2
    adcp(jj).tau_low(kk)=(adcp(jj).ustar(kk).^2)/10;
    adcp(jj).tau_qsl(kk)=adcp(jj).Cd(kk)*(adcp(jj).ubar(kk).^2)/10;
    end
end
tau_low=vertcat(adcp.tau_low);
tau_qsl=vertcat(adcp.tau_qsl);
tau_std=std(tau_qsl,0,2);
tau_med=nanmedian(tau_qsl,2);



for jj=1:3
    figure(1),subplot(211)
    plot(MeasuredTime,tau_qsl(:,jj)),hold on
    figure(2)
    subplot(1,3,jj)
    scatter(tau_low(:,jj),tau_qsl(:,jj)),xlim([0 5])
    refline(1,0)
end
figure(1),datetick('x','HHMM','keeplimits')
subplot(212)
errorbar(MeasuredTime,tau_med,tau_std),hold on
datetick('x','HHMM','keeplimits')