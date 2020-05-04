
%%

% Robin Banner, July 2016, adapted by HEG 7/23/2019
% This code finds the average depth, River Velocity, and Tidal Velocity at 
% a certain trasnect/station and then calculates the freshwater Froude Number and 
% Mixing Number [based on Geyer and MacCready (2014)] and plots the points in 
% estuarine parameter space.

clear all%,close all,clc

So = 27; % Ocean Salinity
Cd=0.0031; %From Sternberg 1968

% F={'YR_Jan19_Neap';'YR_Sept17_Neap';'BR_Mar18_Neap';'BR_Mar18_Spring';'BR_Sept17_Neap';'BR_Sept17_Spring'};
F={'YR_Jan19_Neap';'BR_Mar18_Neap';'BR_Mar18_Spring';'PR_Sept17_Spring'};

% H is determined in advance based on Fatho data
H=[17,16,16,26];

for kk=1:length(F)
    load([F{kk},'.mat'])
    load([F{kk},'_FluxDecomp2.mat'])
    Areas=[];
    vel=[];

% % Mean Thalweg Depth (IGNORE)
% H=vertcat(adcp.interpdepths);
% H=mean(max(H,[],2,'omitnan'))


% River Velocity (Ur = Qr/Area)~=u0 in flux decomp
Ur(kk)=fluxdecomp.U0;
if kk==2
    Ur=abs(Ur)*0.6;
elseif kk==4 % Convert Pathein HF to LF using the ratio of HF/LF in the Yangon
    Ur(kk)=Ur(kk)/10;
end

% Tidal Velocity (amplitude of depth averaged velocity)
for jj=1:length(adcp)
    vel(jj,:)=nanmean(adcp(jj).alongComplete./100);
end
Ut(kk) = (max(max(vel))-min(min(vel)))/2;
if kk==1
    Ut(kk)=Ut(kk)*1.05;
end

% Tidal Excursion (tidal velocity amplitude * tidal period / pi
Ex(kk)=(Ut(kk)*(12.41*60*60))/pi/1000; % m to km

[Frf(kk),M(kk),No(kk)] = eps_calc(H(kk),Ur(kk),Ut(kk),So,Cd);

figure(10)
loglog(M(kk),Frf(kk),'*'),hold on
end
xlim([0.15 3]),ylim([10e-5 1])
ax=gca;ax.XTick=[0.2 0.5 1 2];
legend({'Yangon, Low, Neap';'Bogale, Low, Neap';'Bogale, Low, Spring';...
    'Pathein, Low, Spring'})

% export all the estuary parameters to csv:
A = [H;Ur;Ut;Ex;Frf;M;No];
% writematrix(A,'LF_EstuaryParameters')

%%
function [Frf,M,No] = eps_calc(H,Ur,Ut,So,Cd)
% Constants shared by all times
w = (2*pi)/(12.41*60*60); % [s-1] frequency of semidiurnal tidal cycle
B = 0.00077; % Constant
g = 9.80665; % [m/s] gravity

No = ((B*g*So)/H)^(1/2); % Bouyancy Frequency for max top bottom salinity variation
Frf = Ur/((B*g*So*H)^(1/2)); % Freshwater Froude Number
M = ((Cd*(Ut^2))/(w*No*(H^2)))^(1/2); % Mixing Number

end


