clear all,close all,clc
% sept17 neap covers the full 12.4 hr tidal cycle, no need to add on
names = {'YR_Sept17_Neap';'BR_Sept17_Neap'};
kk=1;
load([names{kk},'.mat'])
load([names{kk},'_SSC.mat'])
load('YangonRiverInstruments.mat')

% get the adcp transect indicies from the ssc transects
sscidx=vertcat(ssc.transect);
% Make a 100x1000 cross section for the sigma format
sigmadepth = linspace(0,1,50).';
L=length(sigmadepth);

for jj=1:length(adcp)
    % find the beginning and end of the transect
    % fix interp depth
    adcp(jj).interpdepths(...
        adcp(jj).interpdepths<0 | adcp(jj).interpdepths>25)=NaN;
    %adcp(jj).interpdepths(50:800)=fillmissing(adcp(jj).interpdepths(50:800),'nearest');
    
    % if there is an ssc matching the current transect, put it in,
    % otherwise, put in the next matching ssc transect
    idx=find(sscidx==jj);
     if isempty(idx)==1
        diff=abs(jj-sscidx);
        [~,idx]=min(diff);
    end
    adcp(jj).Density=nanmean(ssc(idx).Density_avg);
    adcp(jj).ssc=ssc(idx).sigmaAlongComplete./1000; % ssc mg/L in kg/m3 (from mg/l)
    
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
                adcp(jj).sigmaDepth(:,nn),adcp(jj).alongComplete(:,nn),sigmadepth)./100;
        end        %convert from cm/s to m/s
    end
    % fill the sigma velocity to the surface
    adcp(jj).sigmaAlongComplete=fillmissing(adcp(jj).sigmaAlongComplete,'nearest',1);
    % fill missing columns of data
    adcp(jj).sigmaAlongComplete(:,~isnan(adcp(jj).interpdepths))=...
        fillmissing(adcp(jj).sigmaAlongComplete(:,~isnan(adcp(jj).interpdepths))...
        ,'nearest',2);
    
    % downsample the sigma crosssection to reduce noise and enable time
    % extrapolation for incomplete tidal surveys
    xbins=100;ybins=10;
    adcp(jj).sigmaAlongComplete=BlockMean_mean(adcp(jj).sigmaAlongComplete,ybins,xbins);
    adcp(jj).ssc=BlockMean_mean(adcp(jj).ssc,ybins,xbins);
    adcp(jj).interpdepths=BlockMean_mean(adcp(jj).interpdepths,1,xbins);
    
    % calculate total cross-section area for this transect
    A(jj)=xbins*sum(adcp(jj).interpdepths,'omitnan');
    
    % calculate simga bin size (h/z *  m width)
    dA(jj,:)=xbins*(ybins*(adcp(jj).interpdepths/L)); % m*m=m2
     
    % L2006 methods section calculations Fs=Sum(u*c*Abin)
    Fs(jj)=nansum(nansum(...
        adcp(jj).sigmaAlongComplete.*adcp(jj).ssc.*...
        dA(jj,:)))./1000; %m/s * kg/m3 * m2 = kg/s/1000= t/s 
       
    % use Qf and S to calc the u0 and S0 by taking tidally
    % interpolated avg and dividing by avg tidal area
    adcp(jj).Qf=adcp(jj).sigmaAlongComplete.*dA(jj,:);% cm/s to m3/s
    Qf(jj)=nansum(nansum(adcp(jj).Qf));%"integration"
    adcp(jj).S=adcp(jj).ssc.*dA(jj,:);%kg/m3 * m2 = kg/m
    S(jj)=nansum(nansum(adcp(jj).S)); % "integration"
    
    % calculate a time for the transects:
    MeasuredTime(jj)=adcp(jj).time(1);
end

figure;
subplot(311)
plot(YangonRiver.datenum,YangonRiver.Depth,'k')
xlim([MeasuredTime(1)-1/48 MeasuredTime(end)+1/48])
datetick('x','HH:MM','keeplimits')
ylabel('water level')

subplot(3,1,2:3)
yyaxis left
plot(MeasuredTime,Qf,'bo'),hold on
yyaxis right
plot(MeasuredTime,Fs,'ko'),hold on

% interpolate to get the full tidal cycle

InterpTimes=days(...
   [min(MeasuredTime):minutes(1):max(MeasuredTime)]);
vq=MeasuredTime;

A =interp1(vq,A,InterpTimes,'linear');
A0=nanmean(A); %m2

dA = interp1(vq,dA,InterpTimes,'linear');
dA0=nanmean(dA); %m2

Qf=interp1(vq,Qf,InterpTimes,'linear');
U0=nanmean(Qf)/A0; %m/s

S=interp1(vq,S,InterpTimes,'linear');
S0=nanmean(S)/A0; %kg/m3

Fs =interp1(vq,Fs,InterpTimes,'linear');

Di0=nanmean(Qf); % m3/s


Fs0=nanmean(Fs);%t/s
Fr=U0*S0*A0/1000; %t/s


yyaxis left
plot(InterpTimes,Qf,'b-'),hold on
ax=gca; ax.YColor='b';
r=refline(0,0); r.Color='k';
ylim([-0.65e4 2e4]),ylabel('Discharge, (m^3/s)')
yyaxis right
plot(InterpTimes,Fs,'k-'),hold on
ax=gca; ax.YColor='k';
%r=refline(0,0); r.Color='r';
ylim([-2.5 8]),ylabel('Sediment Flux, (t/s)')
xlim([MeasuredTime(1)-1/48 MeasuredTime(end)+1/48])
datetick('x','HH:MM','keeplimits')
