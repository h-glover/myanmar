% calculate Rouse parameter based on profiles of SSC
clear all,close all,clc
fname='YR_Sept17_Neap';
load([fname,'_SSC.mat'],'ssc')

for jj=1:length(ssc)
%     figure;
    cols=find(~isnan(ssc(jj).Depth(1,:)));
    for kk=cols
        % calculate height above bed, assuming max depth reached is the bed
        Z=100*(max(ssc(jj).Depth(:,kk))-ssc(jj).Depth(:,kk))+30;
        
        % find the location of the maximum in SSC for Ca
        maxSSC=find(ssc(jj).SSC(:,kk)==max(ssc(jj).SSC(:,kk)));
        if maxSSC==length(ssc(jj).SSC(:,1))
            Ca(jj,kk)=nanmean(ssc(jj).SSC(maxSSC-2:maxSSC,kk));
        else
            Ca(jj,kk)=nanmean(ssc(jj).SSC(maxSSC-1:maxSSC+2,kk));
        end
        Za=Z(maxSSC);
        
        % Conc profile is the SSC/maxSSC; only use the bottom 5m
        C=ssc(jj).SSC(Z<400,kk)./Ca(jj,kk);
        Z=Z(Z<400)./Za;
        
        %remove nans for the fit calculation and sort into descending order
        Z(isnan(C))=[];
        C(isnan(C))=[];
        [Z,srt]=sort(Z,'ascend');
        C=C(srt);
        
        % fit an exponential curve to the data C = Z^P
        myfittype = fittype('x^a','dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a'});
       
        %w=fliplr(exp(linspace(1,100,length(C))));
        w=ones(1,length(C));
        w([1,15,end])=500;
        
        options=fitoptions(myfittype);
        options.Robust='LAR';
        options.Upper=0;
        options.StartPoint=0.1;
        options.MaxIter=1000;
        options.MaxFunEvals=1000;
        options.TolX=10^(-10);
        options.TolFun=10^(-10);
        options.Weights=w;
        f=fit(Z,C,myfittype,options);
        ssc(jj).P(kk)=-f.a;
%         subplot(1,3,kk)
%         plot(f,Z,C),hold on
        
               
    end
    ssc(jj).P(ssc(jj).P<0.001)=NaN;

end

save([fname,'_SSC.mat'],'ssc')

%% Calculate approx settling velocity from P and ustar
clear all,close all,clc
fname='BR_Sept17_Spring';
F=dir('*_ShearStress.mat');
for ff=6
    fname=F(ff).name(1:end-16);
load([fname,'_SSC.mat'],'ssc')
load([fname,'_ShearStress.mat'])

time_ustar=avgs(1).time;
ustar=cat(3,avgs.ustar);
ustar_std=std(ustar,0,3,'omitnan');
ustar=nanmean(ustar,3);

sscidx=vertcat(ssc.transect);
for jj=1:length(ssc)
    if length(ssc(jj).P)<3
        ssc(jj).P(end+1)=NaN;
    end
end


for jj=1:length(time_ustar)
    idx=find(sscidx==jj);
     if isempty(idx)==1
        diff=abs(jj-sscidx);
        [~,idx]=min(diff);
     end
    P(jj,:)=ssc(idx).P;
end

for jj=1:3
    Ws(:,:,jj)=0.41.*P.*ustar(:,jj);
    Ws(Ws>0.5)=NaN;
end
Ws=nanmedian(Ws,3);
figure;
subplot(1,3,1:2)
plot(time_ustar,Ws)
title(F(ff).name)
datetick('x','HHMM','keeplimits')
ylabel('cm/s settling velocity')
subplot(1,3,3)
histogram(Ws,[0:0.01:0.5])
nanmedian(Ws)

% clearvars -except ff F
end

%% Stokes settling 
