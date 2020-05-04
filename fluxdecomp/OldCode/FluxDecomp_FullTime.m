clear all,clc
% sept17 neap covers the full 12.4 hr tidal cycle, no need to add on
names = {'YR_Sept17_Neap';'BR_Sept17_Neap'};

for kk=1:length(names)
load([names{kk},'.mat'])
load([names{kk},'_SSC.mat'])

% get the adcp transect indicies from the ssc transects
sscidx=vertcat(ssc.transect);
% Make a 100x1000 cross section for the sigma format
sigmadepth = linspace(0,1,100).';
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

    % calculate total cross-section area for this transect
    A(jj)=sum(adcp(jj).interpdepths,'omitnan');
    
    % calculate simga bin size (h/z * 1 m width)
    dA(jj,:)=adcp(jj).interpdepths/L; % m*m=m2
     
    % L2006 methods section calculations Fs=Sum(u*c*Abin)
    Fs(jj)=nansum(nansum(...
        adcp(jj).sigmaAlongComplete.*adcp(jj).ssc.*...
        dA(jj,:)))./1000; %m/s * kg/m3 * m2 = kg/s/1000= t/s 
       
    % calculate total discharge (Fs=Sum(U*Abin))
    Di(jj)=nansum(nansum(...
        adcp(jj).sigmaAlongComplete.*dA(jj,:))); %m/s * m2 = m3/s
    
    % use Qf and S to calc the u0 and S0 by taking tidally
    % interpolated avg and dividing by avg tidal area
    adcp(jj).Qf=adcp(jj).sigmaAlongComplete.*dA(jj,:);% cm/s to m3/s
    Qf(jj)=nansum(nansum(adcp(jj).Qf));%"integration"
    adcp(jj).S=adcp(jj).ssc.*dA(jj,:);%kg/m3 * m2 = kg/m
    S(jj)=nansum(nansum(adcp(jj).S)); % "integration"
    
    % calculate a time for the transects:
    MeasuredTime(jj)=adcp(jj).time(1);
end
fluxdecomp.MeasuredTime=MeasuredTime;
fluxdecomp.Meas.Qf=Qf;
fluxdecomp.Meas.Fs=Fs;
%% interpolate to get the full tidal cycle


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

Di = interp1(vq,Di,InterpTimes,'linear');
Di0=nanmean(Di) % m3/s

names{kk}
Fs0=nanmean(Fs)%t/s
Fr=U0*S0*A0/1000 %t/s

%% calc u1 and s1 by temporally avging then integrating
U1=cat(3,adcp.Qf);
S1=cat(3,adcp.S);

% create meshgrids with the correct times
[X,Y,Z]=meshgrid(1:nn,1:L,vq);
[Xq,Yq,Zq]=meshgrid(1:nn,1:L,InterpTimes);
U1 =interp3(X,Y,Z,U1,Xq,Yq,Zq);
U1=nanmean(U1,3)./dA0 - U0;% m3/s / m2 = m/s

S1 =interp3(X,Y,Z,S1,Xq,Yq,Zq);
S1=nanmean(S1,3)./dA0 - S0;% kg/m / m2 = kg/m3
Fe=nansum(nansum(U1.*S1.*dA0))/1000 %sum(m/s * kg/m3)*m2 =kg/s/1000=t/s


%% u2 and s2
U2=cat(3,adcp.sigmaAlongComplete);% m/s
S2=cat(3,adcp.ssc);% kg/m3

U2 =interp3(X,Y,Z,U2,Xq,Yq,Zq);
U2=U2 - U0 - U1;% m/s

S2 =interp3(X,Y,Z,S2,Xq,Yq,Zq);
S2=S2 - S0 - S1;% kg/m3
dA=reshape(dA',[1,nn,length(A)]);% reshape dA to be multiplied by S2 and U2

Ft=nanmean(nansum(nansum(U2.*S2.*dA,1),2))./1000 % kg/m3*m/s*m2 = kg/s/1000=t/s

%%
fluxdecomp.time=InterpTimes;
fluxdecomp.A=A;
fluxdecomp.A0=A0;
fluxdecomp.dA=dA;
fluxdecomp.dA0=dA0;
fluxdecomp.Fs=Fs;
fluxdecomp.Fs0=Fs0;
fluxdecomp.Di=Di;
fluxdecomp.Di0=Di0;
fluxdecomp.Qf=Qf;
fluxdecomp.U0=U0;
fluxdecomp.S=S;
fluxdecomp.S0=S0;
%fluxdecomp.U1=U1;
%fluxdecomp.S1=S1;
fluxdecomp.Ft=Ft;
fluxdecomp.Fr=Fr;
fluxdecomp.Fe=Fe;

save([names{kk},'_FluxDecomp2'],'fluxdecomp')
clearvars -except names
end