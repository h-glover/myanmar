function [skew_elev_move, skew_vel_move, skew_acc_move,tvec] = skew_asym3(slope,vel,acc,time)

% RLM March 2018
% edited by HEG, 5/18/2020

% slope = d[water level]/dt in m/hr
% vel = water velocity in m/s
% ssc = SSC in g/L
% time = matlab datetime vector

% Skewness
    % I will quantify asymmetry as the normalized sample skewness of the tidal
    % elevation (skew_elevation) and velocity (skew_velocity)time derivative based on the methods of Nidzieko (2010)

    % Skewness>0 : Flood dhdt and velocity duration are shorter
    % Skewness<0 : Ebb dhdt and velocity duration are shorter

    % Main Equation
        % skewness = u3/(sigma^3)
            % u3 = third sample moment about the mean
            % sigma = standard deviation, ie the square root of the second sample moment about the mean

    % Important variables
        % T = number of observations from time=0 to time=T
        % dddt = d(water depth)/d(time)
        % daspd = d(depth averaged velocity)/d(time)
        % moment = sample moment about the mean

% Calculate skewness with a moving window (skew_move)

% calc the time step in hours:
dt=24*(time(2)-time(1));
% moving mean window length for a single tide:
T_move=floor(12.5/dt);

% calculate lenght of time series and make new index vector for sampling at
% 1 : moving mean interval : length of series
L = length(slope);
ind=1:T_move:L; 
tvec = time(ind(1:end-1));
%% calculate asymmetry from the water slope: dh/dt
for n = 1:length(ind)-1
    moment_move = slope(ind(n):ind(n+1)-1);
    u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
    sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
    skew_elev_move(n) = u3_move./(sigma3_move);
end
% % Add last point with backward looking window
% moment_move = slope(L-T_move:L);
% u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
% sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
% skew_elev_move(end+1) = u3_move./(sigma3_move);


%%  Calculate velocity skewness 
for n = 1:length(ind)-1
    moment_move = vel(ind(n):ind(n+1)-1);
    u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
    sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
    skew_vel_move(n) = u3_move./(sigma3_move);
end
% % Add last point with backward looking window
% moment_move = vel(L-T_move:L);
% u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
% sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
% skew_vel_move(end+1) = u3_move./(sigma3_move);



%% --------------------------------------------------------------------------------------
% Calculate skewness with a moving window (skew_move)
% T_move = 1488; % 12.4 hour moving window for 1/60 Hz data

for n = 1:length(ind)-1
    moment_move = acc(ind(n):ind(n+1)-1);
    u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
    sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
    skew_acc_move(n) = u3_move./(sigma3_move);
end
% % Add last point with backward looking window
% moment_move = flux(L-T_move:L);
% u3_move = (1/(T_move-1))*sum(moment_move.^3,'omitnan');
% sigma3_move = ((1/(T_move-1))*sum(moment_move.^2,'omitnan')).^(3/2);
% skew_flux_move(end+1) = u3_move./(sigma3_move);


end
