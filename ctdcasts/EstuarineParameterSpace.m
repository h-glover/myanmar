
%%

% Robin Banner, July 2016, adapted by HEG 7/23/2019
% This code finds the average depth, River Velocity, and Tidal Velocity at 
% a certain trasnect/station and then calculates the freshwater Froude Number and 
% Mixing Number [based on Geyer and MacCready (2014)] and plots the points in 
% estuarine parameter space.

clear all%,close all,clc
% (1)Low Flow Middling

So = 27; % Ocean Salinity
Cd=0.0031; %From Sternberg 1968

F={'YR_Jan19_Neap';'BR_Mar18_Neap';'BR_Mar18_Spring'};

for kk=1:length(F)
    load([F{kk},'.mat'])
    load([F{kk},'_discharge_interp.mat'])
    Areas=[];
    vel=[];
    meandepth_pass=[];
    % thalweg is different in Yangon (first 2) than Bogale
    if kk==1
        thlwg = 175:325;
    elseif kk==2
        InterpRes=150;
        thlwg = 300:600;
    else
        thlwg = 300:600;
    end
    
% Mean Thalweg Depth
for jj = 1:length(adcp)
    adcp(jj).interpdepths(adcp(jj).interpdepths<0)=NaN;
    meandepth_pass(jj) = nanmean(adcp(jj).interpdepths(thlwg));
end
H = nanmean(meandepth_pass);

% River Velocity (Ur = Qr/Area), Qr = net flow over tidal cycle
for jj=1:length(adcp)
    adcp(jj).interpdepths(isnan(adcp(jj).interpdepths))=0;
    Areas(jj)=trapz(adcp(jj).interpdepths);
end
Ur = InterpRes./nanmean(Areas);

% Tidal Velocity (amplitude of depth averaged velocity)
for jj=1:length(adcp)
    %vel(jj)=nanmean(nanmean(adcp(jj).alongComplete./100,2));
    vel(jj)=nanmean(nanmean(abs(adcp(jj).alongComplete./100),2));
    %vel(jj)=nanmean(nanmean(adcp(jj).alongComplete(:,thlwg)./100,2));
end
Ut = (max(vel) - min(vel))/2


[Frf(kk),M(kk),No(kk)] = eps_calc(H,Ur,Ut,So,Cd);
%Flux(kk)=Ur*(12.41*60*60)/H;
Flux(kk)=pi*Ur/Ut;

figure(1)
loglog(M(kk),Frf(kk),'b*'),hold on
end
 
xlim([0.15 3]),ylim([10e-5 1])
ax=gca;ax.XTick=[0.2 0.5 1 2];
%legend({'Yangon, Low, Neap';'Bogale, Low, Neap';'Bogale, Low, Spring'})
%%
clear
H=70;
So = 27; % Ocean Salinity
Cd=0.0031;
w = (2*pi)/(12.41*60*60); % [s-1] frequency of semidiurnal tidal cycle
B = 0.00077; % Constant
g = 9.80665; % [m/s] gravity
No = ((B*g*So)/H)^(1/2); % Bouyancy Frequency for max top bottom salinity variation
Ut = 0.8;
M = ((Cd*(Ut^2))/(w*No*(H^2)))^(1/2) 

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


