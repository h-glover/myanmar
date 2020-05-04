clear all,clc,close all

names = {'BR_Sept17_Neap'};

ff=1;
load([names{ff},'.mat'])
load([names{ff},'_SSC.mat'])

% get the adcp transect indicies from the ssc transects
sscidx=vertcat(ssc.transect);
% Make a 100x1000 cross section for the sigma format
sigmadepth = linspace(0,1,50).';
L=length(sigmadepth);

for jj=1:length(adcp)
    % fix interp depth
    adcp(jj).interpdepths(...
        adcp(jj).interpdepths<0 | adcp(jj).interpdepths>25)=NaN;
        
    % if there is an ssc matching the current transect, put it in,
    % otherwise, put in the next matching ssc transect
    idx=find(sscidx==jj);
     if isempty(idx)==1
        diff=abs(jj-sscidx);
        [~,idx]=min(diff);
    end
    adcp(jj).Density=nanmean(ssc(idx).Density_avg);
    adcp(jj).ssc=ssc(idx).sigmaAlongComplete./1000; % ssc mg/L to kg/m3 
    
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
    adcp(jj).sigmaAlongComplete=fillmissing(adcp(jj).sigmaAlongComplete,...
        'linear',2);
    adcp(jj).ssc=BlockMean_mean(adcp(jj).ssc,ybins,xbins);
    adcp(jj).ssc=abs(fillmissing(adcp(jj).ssc,'linear',2));
    adcp(jj).interpdepths=BlockMean_mean(adcp(jj).interpdepths,1,xbins);
    
    % calculate total cross-section area for this transect
    A(jj)=xbins*sum(adcp(jj).interpdepths,'omitnan');
    
    % calculate simga bin size (h/z *  m width)
    dA(jj,:)=xbins*(ybins*(adcp(jj).interpdepths/L)); % m*m=m2
     
    % L2006 methods section calculations Fs=Sum(u*c*Abin)
    Fs(jj)=nansum(nansum(...
        adcp(jj).sigmaAlongComplete.*adcp(jj).ssc.*...
        dA(jj,:)))./1000; %m/s * kg/m3 = kg/ms * m2 = kg/s/1000= t/s 
    
    % use Qf and S to calc the u0 and S0 by taking tidally
    % interpolated avg and dividing by avg tidal area
    adcp(jj).Qf=adcp(jj).sigmaAlongComplete.*dA(jj,:);% cm/s to m3/s
    Qf(jj)=nansum(nansum(adcp(jj).Qf));%"integration"
    adcp(jj).Qf=fillmissing(adcp(jj).Qf,'linear',2);
    adcp(jj).S=adcp(jj).ssc.*dA(jj,:);%mg/l to kg/m
    S(jj)=nansum(nansum(adcp(jj).S)); % "integration"
    adcp(jj).S=fillmissing(adcp(jj).S,'linear',2);
    
    % calculate a time for the transects:
    MeasuredTime(jj)=adcp(jj).time(1);
end

%% Interp to get full tidal cycle
% Make the first measure the last measure of 12.42hr tidal cycle
MeasuredTime=[MeasuredTime,days(min(MeasuredTime)+hours(12.42))];
Qf=[Qf,Qf(1)];
Fs=[Fs,Fs(1)];
% make the time vector for the full interpolated transect day:
InterpTimes=days(...
   [min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.42)]);

% for x-section area, bin area, and sediment conc do a linear interpolation
% to finish the tidal cycle
A =interp1(MeasuredTime,[A,A(1)],InterpTimes,'linear');
dA = interp1(MeasuredTime,[dA;dA(1,:)],InterpTimes,'linear');
S=interp1(MeasuredTime,[S,S(1)],InterpTimes,'linear');

% for Qf and Fs do a fourier fit and evaluate the rmse for error
[fitobj,gof] = fit(MeasuredTime',Qf','fourier3');
Qf_interp=feval(fitobj,InterpTimes);
Di0=nanmean(Qf_interp);
rmse.Di0=gof.rmse;

figure;subplot(211)
plot(MeasuredTime,Qf,'ko'),hold on,plot(InterpTimes,Qf_interp)

[fitobj,gof] = fit(MeasuredTime',Fs','fourier3');
Fs_interp=feval(fitobj,InterpTimes);
rmse.Fs=gof.rmse;
subplot(212)
plot(MeasuredTime,Fs,'ko'),hold on,plot(InterpTimes,Fs_interp)


% calculate flux values
A0=nanmean(A); %m2
dA0=nanmean(dA); %m2
S0=nanmean(S)/A0; %kg/m3
Di0=nanmean(Qf_interp); % m3/s
U0=Di0/A0; %m/s
Fs0=nanmean(Fs_interp)%t/s
Fr=U0*S0*A0/1000 %t/s

%% Calculate Fe using S1 and U1
% Concatenate the adcp transects, 3rd dim=time
S1=cat(3,adcp.S);
S1=cat(3,S1,adcp(1).S);
[mm,nn,~]=size(S1);

% add first transect to end (fake transect after 12.42hrs)
U1=cat(3,adcp.Qf);
U1=cat(3,U1,adcp(1).Qf);

% create a time series just over the actual measurement period to linearly
% interpolate the actual measuremetns to a 1min interval
LinFillTimes=days([MeasuredTime(1):minutes(1):MeasuredTime(end-1)]);

% preallocate interpolation matrix with NaNs
U1_interp=NaN(mm,nn,length(InterpTimes));
S1_interp=NaN(mm,nn,length(InterpTimes));

for jj=1:mm% double loop for lack of time to find right soln
    for kk=1:nn %make vector for fillmissing func then put back into 3D mat
        temp=reshape(S1(jj,kk,:),[1,length(MeasuredTime)]);
        % fit a spline curve to the fake 25hr timeseries
        [fitobj,gof] = fit(MeasuredTime',temp','fourier3');
        yy=feval(fitobj,InterpTimes);
%         figure;subplot(211)
%         plot(MeasuredTime,temp,'ko'),hold on, plot(InterpTimes,yy,'r')
        rmse.S1(jj,kk)=gof.rmse;
        S1_interp(jj,kk,:)=reshape(yy,[1,1,length(InterpTimes)]);
        
        temp=reshape(U1(jj,kk,:),[1,length(MeasuredTime)]);
        % fit a spline curve to the fake 25hr timeseries
        [fitobj,gof] = fit(MeasuredTime',temp','fourier3');
        yy=feval(fitobj,InterpTimes);
%         subplot(212)
%         plot(MeasuredTime,temp,'ko'),hold on, plot(InterpTimes,yy,'r')
        rmse.U1(jj,kk)=gof.rmse;
        U1_interp(jj,kk,:)=reshape(yy,[1,1,length(InterpTimes)]);
    end
end
rmse.S1=rmse.S1./dA0;
S1=nanmean(S1_interp,3)./dA0 - S0;% kg/m / m2 = kg/m3

rmse.U1=rmse.U1./dA0;
U1=nanmean(U1_interp,3)./dA0 - U0;% m3/s / m2 = m/s

Fe=nansum(nansum(U1.*S1.*dA0))/1000 %sum(m/s * kg/m3)*m2 =kg/s/1000=t/s

Fe_err=sqrt(sum(sum((rmse.U1./U1).^2 + (rmse.S1./S1).^2)))/1000;

%% u2 and s2
U2=cat(3,adcp.sigmaAlongComplete);% m/s
S2=cat(3,adcp.ssc);% kg/m3
U2=cat(3,U2,U2(:,:,1));
S2=cat(3,S2,S2(:,:,1));

U2=U2-U0-U1;
S2=S2-S0-S1;

U2_interp=NaN(mm,nn,length(InterpTimes));
S2_interp=NaN(mm,nn,length(InterpTimes));

for jj=1:mm% double loop for lack of time to find right soln
    for kk=1:nn %make vector for fillmissing func then put back into 3D mat
        temp=reshape(U2(jj,kk,:),[1,length(MeasuredTime)]);
        % fit a spline curve to the fake 25hr timeseries
        [fitobj,gof] = fit(MeasuredTime',temp','fourier3');
        y=feval(fitobj,InterpTimes);
        % and fillnan values using fit output, eval error using rmse
        rmse.U2(jj,kk)=gof.rmse;
        U2_interp(jj,kk,:)=reshape(y,[1,1,length(InterpTimes)]);
        
%         figure; subplot(211),plot(MeasuredTime,temp,'ko'),hold on
%         plot(InterpTimes,y,'r')
        
        temp=reshape(S2(jj,kk,:),[1,length(MeasuredTime)]);
        % fit a spline curve to the fake 25hr timeseries
        [fitobj,gof] = fit(MeasuredTime',temp','fourier3');
        y=feval(fitobj,InterpTimes);
        % and fillnan values using fit output, eval error using rmse
        rmse.S2(jj,kk)=gof.rmse;
        S2_interp(jj,kk,:)=reshape(y,[1,1,length(InterpTimes)]);
        
%         subplot(212),plot(MeasuredTime,temp,'ko'),hold on
%         plot(InterpTimes,y,'r')
        
    end
end
dA=reshape(dA',[1,nn,length(A)]);% reshape dA to be multiplied by S2 and U2
Ft=nanmean(nansum(nansum(U2_interp.*S2_interp.*dA,1),2))./1000; % kg/m3*m/s*m2 = kg/s/1000=t/s

Ft_err=sqrt(sum(sum((...
    rmse.U2./nanmean(U2_interp,3)).^2 + (rmse.S2./nanmean(S2_interp,3)...
    ).^2)))/1000;





%% 
load([names{1},'_FluxDecomp2.mat'])

err=[rmse.Fs,rmse.Fs,Fe_err,Ft_err,0];
figure;
errorbar(1:5,[Fs0,Fr,Fe,Ft,Fs0-(Fr+Fe+Ft)],err,'r')

hold on
plot([fluxdecomp.Fs0,fluxdecomp.Fr,fluxdecomp.Fe,fluxdecomp.Ft,...
    fluxdecomp.Fs0-(fluxdecomp.Fr+fluxdecomp.Fe+fluxdecomp.Ft)],'k')

