clear all,clc
% sept17 neap covers the full 12.5 hr tidal cycle, no need to add on
names = {'BR_Mar18_Spring'};

for kk=1:length(names)
load([names{kk},'.mat'])
load([names{kk},'_Sal.mat'])

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
    adcp(jj).Density=nanmean(sal(idx).Density_avg);
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
        dA(jj,:)))./1000; %m/s * kg/m3 = kg/ms * m2 = kg/s/1000= t/s 
    
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

InterpTimes=days(...
   [min(MeasuredTime):minutes(1):min(MeasuredTime)+hours(12.5)]);
vq=[MeasuredTime,days(min(MeasuredTime)+hours(12.5))];

A =interp1(vq,[A,A(1)],InterpTimes,'linear');
A0=nanmean(A); %m2

dA = interp1(vq,[dA;dA(1,:)],InterpTimes,'linear');
dA0=nanmean(dA); %m2

Qf=interp1(vq,[Qf,Qf(1)],InterpTimes,'linear');
U0=nanmean(Qf)/A0; %m/s

S=interp1(vq,[S,S(1)],InterpTimes,'linear');
S0=nanmean(S)/A0; %kg/m3

Fs =interp1(vq,[Fs,Fs(1)],InterpTimes,'linear');

names{kk}
Fs0=nanmean(Fs)%t/s
Fr=U0*S0*A0/1000 %t/s

%% calc u1 and s1 by temporally avging then integrating
U1=cat(3,adcp.Qf);
S1=cat(3,adcp.S);

% create meshgrids with the correct times
[mm,nn,~]=size(U1);
[X,Y,Z]=meshgrid(1:nn,1:mm,vq);
[Xq,Yq,Zq]=meshgrid(1:nn,1:mm,InterpTimes);
U1(:,:,end+1) = U1(:,:,1);%add first sample on end for interping
U1 =interp3(X,Y,Z,U1,Xq,Yq,Zq);
U1=nanmean(U1,3)./dA0 - U0;% m3/s / m2 = m/s

S1(:,:,end+1) = S1(:,:,1);
S1 =interp3(X,Y,Z,S1,Xq,Yq,Zq);
S1=nanmean(S1,3)./dA0 - S0;% kg/m / m2 = kg/m3
Fe=nansum(nansum(U1.*S1.*dA0))/1000 %sum(m/s * kg/m3)*m2 =kg/s/1000=t/s


%% u2 and s2
U2=cat(3,adcp.sigmaAlongComplete);% m/s
S2=cat(3,adcp.sal);% kg/m3

U2(:,:,end+1) = U2(:,:,1);%add first sample on end for interping
U2 =interp3(X,Y,Z,U2,Xq,Yq,Zq);
U2=U2 - U0 - U1;% m/s

S2(:,:,end+1) = S2(:,:,1);
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
fluxdecomp.Qf=Qf;
fluxdecomp.U0=U0;
fluxdecomp.S=S;
fluxdecomp.S0=S0;
%fluxdecomp.U1=U1;
%fluxdecomp.S1=S1;
fluxdecomp.Ft=Ft;
fluxdecomp.Fr=Fr;
fluxdecomp.Fe=Fe;

figure;
plot(Qf)

save([names{kk},'_SaltFluxDecomp2'],'fluxdecomp')
clearvars -except names
end