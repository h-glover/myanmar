clear all,close all,clc
F=dir('*SS*.mat');

for jj=1:length(F)
    load(F(jj).name)
    if jj==2
        D(11:12)=[];
    elseif jj==1
        D(29:31)=[];
    end
    all(jj).d50=vertcat(D.median);
    all(jj).css=vertcat(D.frac);
	for kk=1:length(D)
        all(jj).data(:,kk)=D(kk).data(:,2);
    end
    all(jj).d50m=mean(all(jj).d50);
end


figure;
subplot(312)
semilogx(D(1).data(:,1),all(1).data,'r'),hold on
plot([all(1).d50m all(1).d50m],[0 5],'r')
semilogx(D(1).data(:,1),all(2).data,'k')
plot([all(2).d50m all(2).d50m],[0 5],'k')
ylabel('Bogale'),xlim([0.3 2000])
ax=gca;ax.XTick=[4 63 1000];

subplot(313)
semilogx(D(1).data(:,1),all(3).data,'r'),hold on
plot([all(3).d50m all(3).d50m],[0 5],'r')
semilogx(D(1).data(:,1),all(4).data,'r')
plot([all(4).d50m all(4).d50m],[0 5],'r')
semilogx(D(1).data(:,1),all(5).data,'k')
plot([all(5).d50m all(5).d50m],[0 5],'k')
ylabel('Pathein'),xlim([0.3 2000])
ax=gca;ax.XTick=[4 63 1000];

subplot(311)
semilogx(D(1).data(:,1),all(6).data,'r'),hold on
plot([all(6).d50m all(6).d50m],[0 5],'r')
semilogx(D(1).data(:,1),all(7).data,'r')
plot([all(7).d50m all(7).d50m],[0 5],'r')
semilogx(D(1).data(:,1),all(8).data,'k')
plot([all(8).d50m all(8).d50m],[0 5],'k')
ylabel('Yangon'),xlim([0.3 2000])
ax=gca;ax.XTick=[4 63 1000];
title('red=Low, black=high')

for jj=1:3
    subplot(3,1,jj)
    plot([10 10],[0 5],'r--')
    plot([30 30],[0 5],'k--')
end