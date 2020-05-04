clear all,clc

names = {'BR_Mar18_Neap';'YR_Jan19_Neap'};

for ff=1%:2
load([names{ff},'.mat'])
load([names{ff},'_Sal.mat'])

% get the adcp transect indicies from the sal transects
salidx=vertcat(sal.transect);
% Make a 100x1000 cross section for the sigma format
sigmadepth = linspace(0,1,50).';
L=length(sigmadepth);

for jj=1:length(adcp)
    % fix interp depth
    adcp(jj).interpdepths(...
        adcp(jj).interpdepths<0 | adcp(jj).interpdepths>25)=NaN;
        
    % if there is an sal matching the current transect, put it in,
    % otherwise, put in the next matching sal transect
    idx=find(salidx==jj);
     if isempty(idx)==1
        diff=abs(jj-salidx);
        [~,idx]=min(diff);
    end
    adcp(jj).Density=nanmean(sal(idx).Density_avg); % 
    % Sal (g/kg=PSU)*kg/m3/1000=kg/m3
    adcp(jj).sal=sal(idx).sigmaAlongComplete.*adcp(jj).Density./1000; 
    
    % put in the CTD locations for future reference
    adcp(jj).CTDlocations=sal(idx).dist;
    
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
    adcp(jj).sigmaAlongComplete=BlockMean(adcp(jj).sigmaAlongComplete,ybins,xbins);
    adcp(jj).sal=BlockMean(adcp(jj).sal,ybins,xbins);
    adcp(jj).interpdepths=BlockMean(adcp(jj).interpdepths,1,xbins);
    
    % calculate total cross-section area for this transect
    A(jj)=xbins*sum(adcp(jj).interpdepths,'omitnan');
    
    % calculate simga bin size (h/z *  m width)
    dA(jj,:)=xbins*(ybins*(adcp(jj).interpdepths/L)); % m*m=m2
     
    % L2006 methods section calculations Fs=Sum(u*c*Abin)
    Fs(jj)=nansum(nansum(...
        adcp(jj).sigmaAlongComplete.*adcp(jj).sal.*...
        dA(jj,:)))./1000; %m/s * kg/m3 * m2 = kg/s/1000= t/s 
    
    % use Qf and S to calc the u0 and S0 by taking tidally
    % interpolated avg and dividing by avg tidal area
    adcp(jj).Qf=adcp(jj).sigmaAlongComplete.*dA(jj,:);% cm/s to m3/s
    Qf(jj)=nansum(nansum(adcp(jj).Qf));%"integration"
    adcp(jj).S=adcp(jj).sal.*dA(jj,:);%kg/m3 to kg/m
    S(jj)=nansum(nansum(adcp(jj).S)); % "integration"
    
    % calculate a time for the transects:
    MeasuredTime(jj)=adcp(jj).time(1);

end
%% interpolate to get the full tidal cycle
names{ff}

% for x-section area, bin area do a linear interpolation
% to finish the tidal cycle
InterpTimes=days(...
   [min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.4)]);
vq=[MeasuredTime,days(min(MeasuredTime)+hours(12.4))];

A =interp1(vq,[A,A(1)],InterpTimes,'linear');
dA = interp1(vq,[dA;dA(1,:)],InterpTimes,'linear');

% for U, S and Fs do a spline fit by repeating the data, making an artificial
% 24hr tidal cycle with a gap in the middle

% first linearly interpolate the existing data
InterpTimes_linear=days([min(MeasuredTime):minutes(1):max(MeasuredTime)]);
InterpQf_linear=interp1(MeasuredTime,Qf,InterpTimes_linear,'linear');
InterpFs_linear=interp1(MeasuredTime,Fs,InterpTimes_linear,'linear');
InterpS_linear=interp1(MeasuredTime,S,InterpTimes_linear,'linear');

% the fill gaps vectors will have [interped_data nans interped_data]
filltimes=days([MeasuredTime(1):minutes(1):MeasuredTime(1)+hours(12.4)]);
InterpTimes_fg=[filltimes,days([filltimes(end)+minutes(1):minutes(1):...
    filltimes(end)+minutes(1)+hours(12.4)])];
LL=length(filltimes)-length(InterpTimes_linear);
InterpQf_fg=[InterpQf_linear,NaN([1,LL]),InterpQf_linear,NaN([1,LL])];
InterpFs_fg=[InterpFs_linear,NaN([1,LL]),InterpFs_linear,NaN([1,LL])];
InterpS_fg=[InterpS_linear,NaN([1,LL]),InterpS_linear,NaN([1,LL])];

%U fit a spline curve to the fake 25hr timeseries
[fitobj,gof] = fit(InterpTimes_fg(~isnan(InterpQf_fg))',...
    InterpQf_fg(~isnan(InterpQf_fg))','fourier1');
y=feval(fitobj,InterpTimes_fg);
figure(100)
plot(fitobj,InterpTimes_fg,InterpQf_fg),title('Qf')
% %and fillnan values using fit output, eval error using rmse
InterpQf_fg(isnan(InterpQf_fg))=y(isnan(InterpQf_fg));
rmse.Qf=gof.rmse;

%Fs fit a spline curve to the fake 25hr timeseries
[fitobj,gof] = fit(InterpTimes_fg(~isnan(InterpFs_fg))',...
    InterpFs_fg(~isnan(InterpFs_fg))','pchipinterp');
y=feval(fitobj,InterpTimes_fg);
% figure(101)
% plot(fitobj,InterpTimes_fg,InterpFs_fg),title('Fs')
% and fillnan values using fit output, eval error using rmse
InterpFs_fg(isnan(InterpFs_fg))=y(isnan(InterpFs_fg));
rmse.Fs=gof.rmse;

%S fit a spline curve to the fake 25hr timeseries
[fitobj,gof] = fit(InterpTimes_fg(~isnan(InterpS_fg))',...
    InterpS_fg(~isnan(InterpS_fg))','pchipinterp');
y=feval(fitobj,InterpTimes_fg);
% figure(101)
% plot(fitobj,InterpTimes_fg,InterpS_fg),title('S')
% and fillnan values using fit output, eval error using rmse
InterpS_fg(isnan(InterpS_fg))=y(isnan(InterpS_fg));
rmse.S=gof.rmse;

% fill using spline and trim to 12.5hr timeseries
Qf=InterpQf_fg(1:end-length(filltimes));
Fs=InterpFs_fg(1:end-length(filltimes));
S=InterpS_fg(1:end-length(filltimes));

% calculate flux values
A0=nanmean(A); %m2
dA0=nanmean(dA); %m2
S0=nanmean(S)/A0; %kg/m3
Di0=nanmean(Qf) % m3/s
U0=nanmean(Qf)/A0; %m/s
Fs0=nanmean(Fs)%t/s
Fr=U0*S0*A0/1000 %t/s

%% calc u1 and s1 by temporally avging then integrating
U1=cat(3,adcp.Qf);
S1=cat(3,adcp.S);

% U: create a "fillgap" matrix like above and fill using spline
% first linearly interpolate for measurement period
[mm,nn,~]=size(U1);
[X,Y,Z]=meshgrid(1:nn,1:mm,MeasuredTime);
[Xq,Yq,Zq]=meshgrid(1:nn,1:mm,InterpTimes_linear);
U1 =interp3(X,Y,Z,U1,Xq,Yq,Zq);
% make 25hr fake time series with gap in middle
U1_fg=cat(3,U1,NaN([mm,nn,LL]),U1,NaN([mm,nn,LL]));

% double loop for lack of time to find right soln
for jj=1:mm
    for kk=1:nn %make vector for fillmissing func then put back into 3D mat
        temp=reshape(U1_fg(jj,kk,:),[1,length(InterpTimes_fg)]);
        % fit a spline curve to the fake 25hr timeseries
        [fitobj,gof] = fit(InterpTimes_fg(~isnan(temp))',...
            temp(~isnan(temp))','pchipinterp');
        y=feval(fitobj,InterpTimes_fg);
        % and fillnan values using fit output, eval error using rmse
        temp(isnan(temp))=y(isnan(temp));
        rmse.U1(jj,kk)=gof.rmse;
        U1_fg(jj,kk,:)=reshape(temp,[1,1,length(InterpTimes_fg)]);
    end
end
U1=U1_fg(:,:,1:end-length(filltimes)); %trim off fake half of time series
U1=nanmean(U1,3)./dA0 - U0;% m3/s / m2 = m/s

% S: create meshgrids with the correct times
S1 =interp3(X,Y,Z,S1,Xq,Yq,Zq);
% make 25hr fake time series with gap in middle
S1_fg=cat(3,S1,NaN([mm,nn,LL]),S1,NaN([mm,nn,LL]));
% double loop for lack of time to find right soln
for jj=1:mm
    for kk=1:nn %make vector for fillmissing func then put back into 3D mat
        temp=reshape(S1_fg(jj,kk,:),[1,length(InterpTimes_fg)]);
        % fit a spline curve to the fake 25hr timeseries
        [fitobj,gof] = fit(InterpTimes_fg(~isnan(temp))',...
            temp(~isnan(temp))','pchipinterp');
        y=feval(fitobj,InterpTimes_fg);
        % and fillnan values using fit output, eval error using rmse
        temp(isnan(temp))=y(isnan(temp));
        rmse.S1(jj,kk)=gof.rmse;
        S1_fg(jj,kk,:)=reshape(temp,[1,1,length(InterpTimes_fg)]);
    end
end
S1=S1_fg(:,:,1:end-length(filltimes)); %trim off fake half of time series
S1=nanmean(S1,3)./dA0 - S0;% kg/m / m2 = kg/m3



Fe=nansum(nansum(U1.*S1.*dA0))/1000 %sum(m/s * kg/m3)*m2 =kg/s/1000=t/s
% %%
% testU=cat(3,adcp.Qf);
% testS=cat(3,adcp.S);
% for jj=1:mm
%     for kk=1:nn
%         dataplot1=reshape(testU(jj,kk,:),[1,length(MeasuredTime)]);
%         dataplot2=reshape(U1(jj,kk,:),[1,length(InterpTimes)]);
%         figure(1)
%         plot(MeasuredTime,dataplot1,'*'),hold on,plot(InterpTimes,dataplot2,'-')
%         
%         dataplot1=reshape(testS(jj,kk,:),[1,length(MeasuredTime)]);
%         dataplot2=reshape(S1(jj,kk,:),[1,length(InterpTimes)]);
%         figure(2)
%         plot(MeasuredTime,dataplot1,'*'),hold on,plot(InterpTimes,dataplot2,'-')
%     end
% end


%% u2 and s2
U2=cat(3,adcp.sigmaAlongComplete);% m/s
S2=cat(3,adcp.sal);% kg/m3

% fill for the measured time using linear interpolation
U2 =interp3(X,Y,Z,U2,Xq,Yq,Zq);
% make 25hr fake time series with gap in middle
U2_fg=cat(3,U2,NaN([mm,nn,LL]),U2,NaN([mm,nn,LL]));
% double loop for lack of time to find right soln for fillmissing in 3D
%{later try permute, reshape, fillmissing(window)?}
for jj=1:mm
    for kk=1:nn %make vector for fillmissing func then put back into 3D mat
        temp=reshape(U2_fg(jj,kk,:),[1,length(InterpTimes_fg)]);
        % fit a spline curve to the fake 25hr timeseries
        [fitobj,gof] = fit(InterpTimes_fg(~isnan(temp))',...
            temp(~isnan(temp))','pchipinterp');
        y=feval(fitobj,InterpTimes_fg);
        % and fillnan values using fit output, eval error using rmse
        temp(isnan(temp))=y(isnan(temp));
        rmse.U2(jj,kk)=gof.rmse;
        U2_fg(jj,kk,:)=reshape(temp,[1,1,length(InterpTimes_fg)]);
    end
end
U2=U2_fg(:,:,1:end-length(filltimes)); %trim off fake half of time series
U2=U2 - U0 - U1;% m/s

%
S2 =interp3(X,Y,Z,S2,Xq,Yq,Zq); %linearly interpolate
% make 25hr fake time series with gap in middle
S2_fg=cat(3,S2,NaN([mm,nn,LL]),S2,NaN([mm,nn,LL]));
% double loop for lack of time to find right soln for fillmissing in 3D
%{later try permute, reshape, fillmissing(window)?}
for jj=1:mm
    for kk=1:nn %make vector for fillmissing func then put back into 3D mat
        temp=reshape(S2_fg(jj,kk,:),[1,length(InterpTimes_fg)]);
        % fit a spline curve to the fake 25hr timeseries
        [fitobj,gof] = fit(InterpTimes_fg(~isnan(temp))',...
            temp(~isnan(temp))','pchipinterp');
        y=feval(fitobj,InterpTimes_fg);
        % and fillnan values using fit output, eval error using rmse
        temp(isnan(temp))=y(isnan(temp));
        rmse.S2(jj,kk)=gof.rmse;
        S2_fg(jj,kk,:)=reshape(temp,[1,1,length(InterpTimes_fg)]);
    end
end
S2=S2_fg(:,:,1:end-length(filltimes)); %trim off fake half of time series
S2=S2 - S0 - S1;% kg/m3


dA=reshape(dA',[1,nn,length(A)]);% reshape dA to be multiplied by S2 and U2
Ft=nanmean(nansum(nansum(U2.*S2.*dA,1),2))./1000 % kg/m3*m/s*m2 = kg/s/1000=t/s

% %%
% testU=cat(3,adcp.sigmaAlongComplete);% m/s
% testS=cat(3,adcp.sal);
% for jj=1:mm
%     for kk=1:nn
%         dataplot1=reshape(testU(jj,kk,:),[1,length(MeasuredTime)]);
%         dataplot2=reshape(U2(jj,kk,:),[1,length(InterpTimes)]);
%         figure(10)
%         plot(MeasuredTime,dataplot1,'*'),hold on,plot(InterpTimes,dataplot2,'-')
%         
%         dataplot1=reshape(testS(jj,kk,:),[1,length(MeasuredTime)]);
%         dataplot2=reshape(S2(jj,kk,:),[1,length(InterpTimes)]);
%         figure(20)
%         plot(MeasuredTime,dataplot1,'*'),hold on,plot(InterpTimes,dataplot2,'-')
%     end
% end
%%
fluxdecomp.time=InterpTimes;
fluxdecomp.A=A;
fluxdecomp.A0=A0;
fluxdecomp.dA=dA;
fluxdecomp.dA0=dA0;
fluxdecomp.Fs=Fs;
fluxdecomp.Fs0=Fs0;
fluxdecomp.Qf=Qf;
fluxdecomp.Di0=Di0;
fluxdecomp.U0=U0;
fluxdecomp.S=S;
fluxdecomp.S0=S0;
%fluxdecomp.U1=U1;
%fluxdecomp.S1=S1;
fluxdecomp.Ft=Ft;
fluxdecomp.Fr=Fr;
fluxdecomp.Fe=Fe;

%save([names{ff},'_SaltFluxDecomp2'],'fluxdecomp')
%clearvars -except names
end

%%
% %% interpolate to get the full tidal cycle
% 
% InterpTimes=days(...
%    [min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.4)]);
% vq=[MeasuredTime,days(min(MeasuredTime)+hours(12.4))];
% 
% A =interp1(vq,[A,A(1)],InterpTimes,'linear');
% A0=nanmean(A); %m2
% 
% dA = interp1(vq,[dA;dA(1,:)],InterpTimes,'linear');
% dA0=nanmean(dA); %m2
% 
% % Qf=interp1(vq,[Qf,Qf(1)],InterpTimes,'linear');
% % U0=nanmean(Qf)/A0; %m/s
% % 
% % S=interp1(vq,[S,S(1)],InterpTimes,'linear');
% % S0=nanmean(S)/A0; %kg/m3
% % 
% Fs =interp1(vq,[Fs,Fs(1)],InterpTimes,'linear');
% 
% names{kk}
% Fs0=nanmean(Fs)%t/s
% % Fr=U0*S0*A0/1000 %t/s
% % 
% 
% InterpTimes_linear=days([min(MeasuredTime):minutes(1):max(MeasuredTime)]);
% InterpQf_linear=interp1(MeasuredTime,Qf,InterpTimes_linear,'linear');
% InterpS_linear=interp1(MeasuredTime,S,InterpTimes_linear,'linear');
% %InterpFs_linear=interp1(MeasuredTime,Fs,InterpTimes_linear,'linear');
% 
% % the fill gaps will have interptimes nans interptimes
% filltimes=days([MeasuredTime(1):minutes(1):MeasuredTime(1)+hours(12.4)]);
% InterpTimes_fg=[filltimes,days([filltimes(end)+minutes(1):minutes(1):...
%     filltimes(end)+minutes(1)+hours(12.4)])];
% LL=length(filltimes)-length(InterpTimes_linear);
% InterpTimes=InterpTimes_fg(1:end-length(filltimes));
% 
% InterpQf_fg=[InterpQf_linear,NaN([1,LL]),InterpQf_linear,NaN([1,LL])];
% InterpQf_fg=fillmissing(InterpQf_fg,'linear');
% InterpQf=InterpQf_fg(1:end-length(filltimes));
% InterpS_fg=[InterpS_linear,NaN([1,LL]),InterpS_linear,NaN([1,LL])];
% InterpS_fg=fillmissing(InterpS_fg,'linear');
% InterpS=InterpS_fg(1:end-length(filltimes));
% % InterpFs_fg=[InterpFs_linear,NaN([1,LL]),InterpFs_linear,NaN([1,LL])];
% % InterpFs_fg=fillmissing(InterpFs_fg,'spline');
% % InterpFs=InterpFs_fg(1:end-length(filltimes));
% 
% 
% Di=nanmean(InterpQf)
% U0=Di/A0;
% 
% S0=nanmean(InterpQf)/A0;
% %Fs0=nanmean(Fs)
% Fr=U0*S0*A0/1000 %t/s
% 
% figure;
% yyaxis left
% plot(InterpTimes,InterpQf),refline(0,Di)
% hold on
% plot(MeasuredTime,Qf,'*')
% yyaxis right
% plot(InterpTimes,Fs),refline(0,Fs0),ylim([-40 40])
% 
% %% calc u1 and s1 by temporally avging then integrating
% U1=cat(3,adcp.Qf);
% S1=cat(3,adcp.S);
% 
% % create meshgrids with the correct times
% [X,Y,Z]=meshgrid(1:nn,1:L,vq);
% [Xq,Yq,Zq]=meshgrid(1:nn,1:L,InterpTimes);
% U1(:,:,end+1) = U1(:,:,1);%add first sample on end for interping
% U1 =interp3(X,Y,Z,U1,Xq,Yq,Zq);
% U1=nanmean(U1,3)./dA0 - U0;% m3/s / m2 = m/s
% 
% S1(:,:,end+1) = S1(:,:,1);
% S1 =interp3(X,Y,Z,S1,Xq,Yq,Zq);
% S1=nanmean(S1,3)./dA0 - S0;% kg/m / m2 = kg/m3
% Fe=nansum(nansum(U1.*S1.*dA0))/1000 %sum(m/s * kg/m3)*m2 =kg/s/1000=t/s
% 
% 
% %% u2 and s2
% U2=cat(3,adcp.sigmaAlongComplete);% m/s
% S2=cat(3,adcp.sal);% kg/m3
% 
% U2(:,:,end+1) = U2(:,:,1);%add first sample on end for interping
% U2 =interp3(X,Y,Z,U2,Xq,Yq,Zq);
% U2=U2 - U0 - U1;% m/s
% 
% S2(:,:,end+1) = S2(:,:,1);
% S2 =interp3(X,Y,Z,S2,Xq,Yq,Zq);
% S2=S2 - S0 - S1;% kg/m3
% dA=reshape(dA',[1,nn,length(A)]);% reshape dA to be multiplied by S2 and U2
% 
% Ft=nanmean(nansum(nansum(U2.*S2.*dA,1),2))./1000 % kg/m3*m/s*m2 = kg/s/1000=t/s
% end
% 
% %% bar chart
% diff=Fs0- (Fe+Fr+Ft);
% M=[Fs0;Fr;Fe;Ft;diff];
% 
% figure;
% subplot(121)
% bar(M)
% load([names{kk},'_SaltFluxDecomp.mat'])
% diff=Fs0- (Fe+Fr+Ft);
% M_lin=[fluxdecomp.Fs0;fluxdecomp.Fr;fluxdecomp.Fe;fluxdecomp.Ft;diff];
% 
% subplot(122)
% bar(M_lin)
% %%
% load('YangonWL_Prediction.mat')
% p_time=p_time+(6.5/24);
% wl=interp1(p_time,YOUT2,InterpTimes_linear,'linear');
% wl2=interp1(p_time,YOUT2,InterpTimes,'linear');
% 
% figure(40)
% scatter(wl,InterpQf_linear),hold on
% scatter(wl2,InterpQf)

