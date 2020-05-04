% plot Pathein water level
clear all,close all,clc
F=dir('Gue*.csv');

for jj=1:length(F)
    T=[];
    fid=fopen(F(jj).name);
    T=textscan(fid,'%*f %s %f %*f','Delimiter',',','HeaderLines',2);
    fclose(fid);
    pr(jj).wl_time=datenum(T{1},'mm/dd/yy HH:MM:SS AM');
    pr(jj).wl=T{2}/100;
    
end

tvec=vertcat(pr.wl_time);
wl=vertcat(pr.wl);
wl(wl<11)=NaN;

t_interp=tvec(1):1/(24*60):tvec(end);
wl_interp=NaN([1,length(t_interp)]);

[~,id,ii]=intersect(tvec,t_interp);
wl_interp(ii)=wl(id);
wl_interp=fillmissing(wl_interp,'linear');

save('PR_WL','wl_interp','t_interp')
%%
coef=ut_solv(t_interp,wl_interp,[],16,'auto','notrend');
[sl_fit,~]=ut_reconstr(t_interp,coef);
figure;
plot(t_interp,wl_interp,'k'),hold on
plot(t_interp,sl_fit,'r')

