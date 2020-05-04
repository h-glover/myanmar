clear all,%close all,clc

load('AyeMar18_CTD_all.mat')
ctd=CombinedProfiles;clear CombinedProfiles
load('CombinedCTDSolo.mat')
solo.turb(solo.turb<0)=0;
for jj=1:length(ctd)
    timediff=abs(ctd(jj).time(1)-solo.time);
    ctdtimes(jj)=ctd(jj).time(1);
    [timediff,idx]=sort(timediff,'ascend');
    if timediff(1)<(3/(60*24))
        idx=idx(1)-1:idx(1)+1;
        ctd(jj).solo=nanmean(solo.turb(idx)).*0.726;
    else
        ctd(jj).solo=NaN;
    end
    
    if isempty(ctd(jj).SSC)==1
        ctd(jj).SSC=NaN;
    end
end
ssc=vertcat(ctd.SSC);
soloturb=vertcat(ctd.solo);
soloturb=soloturb(ssc<1300);ssc=ssc(ssc<1300);
[ssc,sortssc]=sort(ssc);
soloturb=soloturb(sortssc);

% calculate a best fit line 
p = polyfit(soloturb,ssc, 1);
yfit1=polyval(p,soloturb);
R1 = corrcoef(yfit1,ssc);
R1= R1(1,2).^2;

% calculate a best fit line with zero intercept
p2= soloturb\ssc;
yfit2 = p2*soloturb;
R2 = corrcoef(ssc,yfit2);
R2= R2(1,2).^2;

figure;
scatter(soloturb,ssc),hold on
plot(soloturb, yfit1, 'k')
plot(soloturb, yfit2,'k--')
t1=['- y=',num2str(p(1)),'x + ',num2str(p(2)),'\newline r^2=',num2str(R1)];
text(700,150,t1)
t2=['-- y=',num2str(p2),'x','\newline r^2=',num2str(R2)];
text(700,370,t2)


%% try calibrating the solo to the OBS then converting...
clear all,close all,clc

load('AyeMar18_CTD_all.mat')
ctd=CombinedProfiles;clear CombinedProfiles
load('CombinedCTDSolo.mat')

for jj=1:length(ctd)
    time(jj)=ctd(jj).time(1);
end
[~,srt]=sort(time,'ascend');
ctd=ctd(srt);
time=vertcat(ctd.time);
time=datevec(time);
time(:,6)=round(time(:,6));
time=datenum(time);

ctdturb=vertcat(ctd.SSCCal);

% time_interp=time(1):datenum(0,0,0,0,0,1):time(end);

% figure;
% plot(time,turb,'r'),hold on
% plot(solo.time,solo.turb,'k')

[~,isolo,ictd]=intersect(solo.time,time);
solo.turb=solo.turb(isolo);
ctdturb=ctdturb(ictd);
solo.turb(isnan(ctdturb))=[];
ctdturb(isnan(ctdturb))=[];

% solo.turb(ctdturb>1500 | isnan(ctdturb))=[];
% ctdturb(ctdturb>1500 | isnan(ctdturb))=[];
% ctdturb(solo.turb>1500)=[];
% solo.turb(solo.turb>1500)=[];


% calculate a best fit line 
p = polyfit(solo.turb,ctdturb, 1)
yfit1=polyval(p,solo.turb);
R1 = corrcoef(yfit1,solo.turb);
R1= R1(1,2).^2;

% calculate a best fit line with zero intercept
p2= solo.turb\ctdturb
yfit2 = p2*solo.turb;
R2 = corrcoef(ctdturb,yfit2);
R2= R2(1,2).^2;

figure;
scatter(solo.turb,ctdturb),hold on
plot(solo.turb, yfit1, 'g')
plot(solo.turb, yfit2,'r')
% t1=['- y=',num2str(p(1)),'x + ',num2str(p(2)),'\newline r^2=',num2str(R1)];
% text(700,150,t1)
% t2=['-- y=',num2str(p2),'x','\newline r^2=',num2str(R2)];
% text(700,370,t2)

% SOLO CALIBRATION IS SSC = solo*0.726
%% Compile Solo data into 1 vector
% F=dir('07*.txt');
% solo.time=[];
% solo.turb=[];
% for jj=1:length(F)
%     temp=[];
%     fid=fopen(F(jj).name);
%     T=textscan(fid,'%s %f','Delimiter',',','HeaderLines',1);
%     temp=datenum(T{1},'yyyy-mm-dd HH:MM:SS.FFF');
%     solo.time=[solo.time;temp];
%     solo.turb=[solo.turb;T{2}];
% end
% 
% [solo.time,idx]=sort(solo.time,'ascend');
% solo.turb=solo.turb(idx);
% 
% save('CombinedCTDSolo','solo')