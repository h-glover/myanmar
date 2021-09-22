close all
ff=dir('*_FluxDecomp3.mat');

for jj=1:length(ff)
    load(ff(jj).name)
    [m,~]=size(fluxdecomp.Qf);
    if m>1
        u=fluxdecomp.Qf'./fluxdecomp.A;
    else
        u=fluxdecomp.Qf./fluxdecomp.A;
    end
    
    mu3=(1/length(u))*sum((u.^3));
    mu2=(1/length(u))*sum((u.^2));
    
    gamma0=mu3./(mu2.^(2/3))
    
    figure;
    subplot(211)
    plot(u),refline(0,0)
    title(ff(jj).name)
    subplot(212)
    histogram(fluxdecomp.Qf'./fluxdecomp.A,[-1.3:0.1:1.3])
    
end

%% compare all velocity values?

%% compare water level speed
close all,clear all
load('YangonRiverInstruments.mat')
deltaT=round((1/24)./(YangonRiver.datenum(2)-YangonRiver.datenum(1)));
slp=[YangonRiver.Depth(2:end)-YangonRiver.Depth(1:end-1),NaN].*deltaT;
slp(abs(slp)>4)=NaN;
% slp(isnan(slp))=[];
t_cycles=[1:round(12.4*60/10):length(slp)];
for jj=1:length(t_cycles)-1
    tvec(jj)=nanmean(YangonRiver.datenum(t_cycles(jj):t_cycles(jj+1)));
    V=[];
    V=slp(t_cycles(jj):t_cycles(jj+1));
    V=V(~isnan(V));
    mu3(jj)=(1/length(V))*sum((V.^3));
    mu2(jj)=(1/length(V))*sum((V.^2));
end
gamma0=mu3./(mu2.^(2/3));

figure; plot(tvec,gamma0),datetick

% divide into seasons