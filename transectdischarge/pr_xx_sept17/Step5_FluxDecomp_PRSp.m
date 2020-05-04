% combine ssc, adcp, and water level
% extrapolate discharge and get a representative tidal velocity
clear all,close all,clc
load('PR_WL.mat')
load('PR_Sept17_Spring')
load('PR_Sept17_Spring_SSC.mat')

% calculate discharge for each transect:
% get the adcp transect indicies from the ssc transects
sscidx=vertcat(ssc.transect);
% Make a 100x1000 cross section for the sigma format
sigmadepth = linspace(0,1,100).';
L=length(sigmadepth);

for jj=1:length(adcp)
    % fix interp depth
    adcp(jj).interpdepths(...
        adcp(jj).interpdepths<0 | adcp(jj).interpdepths>40)=NaN;
        
    % if there is an ssc matching the current transect, put it in,
    % otherwise, put in the next matching ssc transect
    idx=find(sscidx==jj);
    if jj==4
        idx=3;
    elseif jj==6
        idx=5;
    elseif jj==8
        idx=6;
    end
    adcp(jj).ssc=ssc(idx).sigmaAlongComplete;
   
    % put in the CTD locations for future reference
    adcp(jj).CTDlocations=ssc(idx).dist;
    
    % divide adcp depths by the max depth
    adcp(jj).sigmaDepth=adcp(jj).zComplete./adcp(jj).interpdepths;

    % for each column with depth>0, interpolate to the sigma depth
    % vector (100 linearly spaced pts between 0:1)
    for nn=1:length(adcp(jj).interpdepths)
        if isnan(adcp(jj).interpdepths(nn))==1
            adcp(jj).sigmaAlongComplete(:,nn)=NaN([L,1]);
        else
            adcp(jj).sigmaAlongComplete(:,nn)=interp1(...
                adcp(jj).sigmaDepth(:,nn),adcp(jj).alongComplete(:,nn),sigmadepth);
        end       
    end
    % fill the sigma velocity to the surface
    adcp(jj).sigmaAlongComplete=fillmissing(adcp(jj).sigmaAlongComplete,'nearest',1);

    % calculate total cross-section area for this transect
    A(jj)=sum(adcp(jj).interpdepths,'omitnan');
    
    % calculate simga bin size (h/z * 1 m width)
    dA(jj,:)=adcp(jj).interpdepths/L; % z/h * 1m = m
     
    % L2006 methods section calculations Fs=Sum(u*c*Abin)
    Fs(jj)=nansum(nansum(...
        (adcp(jj).sigmaAlongComplete./100).*(adcp(jj).ssc./1000).*...
        dA(jj,:)))./1000; %m/s * kg/m3 = kg/ms * m2 = kg/s/1000= t/s 
    
    % use Qf and S to calc the u0 and S0 by taking tidally
    % interpolated avg and dividing by avg tidal area
    adcp(jj).Qf=(adcp(jj).sigmaAlongComplete./100).*dA(jj,:);% cm/s to m3/s
    Qf(jj)=nansum(nansum(adcp(jj).Qf));%"integration"
    adcp(jj).S=(adcp(jj).ssc./1000).*dA(jj,:);%mg/l to kg/m
    S(jj)=nansum(nansum(adcp(jj).S)); % "integration"
    
    % calculate a time for the transects:
    MeasuredTime(1,jj)=adcp(jj).time(1);
    

end

MeasuredTime=datevec(MeasuredTime);
MeasuredTime(:,6)=0;
MeasuredTime=datenum(MeasuredTime);
%
slp=wl_interp(2:end)-wl_interp(1:end-1);slp(end+1)=slp(end);

[~,is,ia]=intersect(t_interp,MeasuredTime);

pQf=polyfit(slp(is),Qf,3);
pFs=polyfit(slp(is),Fs,3);
Qffit=polyval(pQf,slp);
Qferr=polyval(pQf,slp(is));
err=sqrt(immse(Qf,Qferr));
Fsfit=polyval(pFs,slp);

figure;
subplot(121)
scatter(slp(is),Qf),hold on
plot(slp,Qffit,'.')
subplot(122)
scatter(slp(is),Fs),hold on
plot(slp,Fsfit,'.')

idx=find(t_interp==MeasuredTime(1));
Qffit=Qffit(idx:idx+(24.8*60));
Fsfit=Fsfit(idx:idx+(24.8*60));
Qffit(Qffit<-10000)=-10000;
Fsfit(Fsfit<-2)=-2;


t_cycle=t_interp(idx:idx+(24.8*60));
figure;
subplot(121)
plot(t_cycle,Qffit),hold on
plot(MeasuredTime,Qf,'*')
subplot(122)
plot(t_cycle,Fsfit),hold on
plot(MeasuredTime,Fs,'*'),ylim([-5.1 5.1])


A0=nanmean(A);
Qf0=nanmean(Qffit)
Fs0=nanmean(Fsfit)
Ur=Qf0./A0
fluxdecomp.time=t_cycle;
fluxdecomp.A=A;
fluxdecomp.A0=A0;
fluxdecomp.dA=dA;
%fluxdecomp.dA0=dA0;
fluxdecomp.Fs=Fs;
fluxdecomp.Fs0=Fs0;
fluxdecomp.Qf=Qf;
fluxdecomp.U0=Ur;


save('PR_Sept17_Spring_FluxDecomp2','fluxdecomp')
