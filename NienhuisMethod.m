% Bogale
clear all
wu = 1000;
L = 17/0.00004;
a = 2.75;
k = 1.1e-4;
b = wu/17;

wm = wu + (b*k*a*L)

%%
clear all,close all,clc

fid=fopen('RawData_short.csv');
rawdat=textscan(fid,'%s %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',2,...
    'Delimiter',',');
fclose(fid);

T=logspace(-4,4,30);
U=[0.1,1,10,100];
for jj=1:length(U)
    ww(jj,:)=T/U(jj)+1;
end

figure;
subplot(121)
scatter(rawdat{5},rawdat{6},'ko'),hold on
scatter(rawdat{5}(2:12),rawdat{6}(2:12),'r*')
scatter(rawdat{5}(1),rawdat{6}(1),'b*')
scatter(rawdat{5}([3,6,9,12]),rawdat{6}([3,6,9,12]),'g*')
xlabel('Obs width'),ylabel('Pred width')
ax=gca; ax.YScale='log';ax.XScale='log';
axis([1 10^5 1 10^5])
refline(1,0)
legend({'All Deltas','Mekong Channels','Mekong','Mekong Distributaries'})

subplot(122)
scatter(rawdat{8},rawdat{10},'ko'),hold on
scatter(rawdat{8}(2:12),rawdat{10}(2:12),'r*')
scatter(rawdat{8}(1),rawdat{10}(1),'b*')
scatter(rawdat{8}([3,6,9,12]),rawdat{10}([3,6,9,12]),'g*')
for jj=1:length(U)
    plot(T,ww(jj,:),'Color',[0.6,0.6,0.6])
end
plot([1 1],[0.5 100],'k--')
ax=gca; ax.YScale='log';ax.XScale='log';
axis([10e-4 10e4 0.5 100])
R=refline(0,1);R.LineStyle='--';R.Color='k';
xlabel('T=Qt/Qr'),ylabel('Wm / Wu')