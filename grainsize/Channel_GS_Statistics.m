clear all,close all,clc

cd C:\GLOVER\output\myanmar\grainsize
load('PR_GrainSize_Mar18.mat')


D([1,3,5,10])=[];
frac = vertcat(D.frac);


p_csv=vertcat(D.lat);
p_csv(:,2)=vertcat(D.long);
p_csv(:,3)=vertcat(D.median);
p_csv(:,5:7)=vertcat(D.frac);

mud = nanmean([p_csv(:,5)+p_csv(:,6)]);

figure;
plot(p_csv(:,2),p_csv(:,1),'o')

% sept 42% mud

%%
csvwrite('GrainSize_Pathein_LowFlow.csv',p_csv)