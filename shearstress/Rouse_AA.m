% This code uses Law of the Wall and Quadratic Stress Law to find U* and Ws

clear all
close all
clc

load('Ased.mat');
load('CTDAdistance.mat');

%% Order by time
for q=1:length(Ased);
    Ased(q).LTmean=mean(Ased(q).DateNumberLT);
end
    timevec=vertcat(Ased(1).LTmean,Ased(2).LTmean,Ased(3).LTmean,Ased(4).LTmean,Ased(5).LTmean,Ased(6).LTmean,Ased(7).LTmean,Ased(8).LTmean,Ased(9).LTmean,Ased(10).LTmean,Ased(11).LTmean,Ased(12).LTmean,Ased(13).LTmean, Ased(14).LTmean,Ased(15).LTmean,Ased(16).LTmean,Ased(17).LTmean,Ased(18).LTmean,Ased(19).LTmean,Ased(20).LTmean,Ased(21).LTmean,Ased(22).LTmean,Ased(23).LTmean,Ased(24).LTmean);
    newtimevec=sort(timevec);
    
for i=1:length(Ased);
    for j=1:length(Ased);
        if newtimevec(i)==timevec(j);
            Asedsort(i,1)=Ased(j);
        end
    end
end
Ased = Asedsort;

%% Average velocities based on sigma coordinates
Station = CTDAdistance;
Stationdepth = 9.7;
sigint = linspace(0,1,50);

for n = 1:length(Ased);
    x = (Ased(n).distance < Station+50).*Ased(n).distance;
    xx = (x > Station-50).*Ased(n).distance;
    xx(isnan(xx))=0;
    I = find (xx ~= 0);
    Ased(n).sigmacrop = ones(Ased(n).numcells,length(Ased(n).sigma));
    for m = I;
        Ased(n).ACvel(isnan(Ased(n).ACvel)) = 0;
        Ased(n).sigma(isnan(Ased(n).sigma)) = 0;
        if Ased(n).ACvel(:,m) == 0;
            continue,
        else
        index = find(Ased(n).ACvel(:,m));
        index = index(end);
        for d = index+1:Ased(n).numcells;
            Ased(n).sigmacrop(d,m) = 0;
        end
        Ased(n).sigmacrop(:,m) = Ased(n).sigmacrop(:,m).* Ased(n).sigma(:,m);
        Ased(n).sigmacrop(Ased(n).sigmacrop(1:index,m)==0,m) = 1;
        Ased(n).velinterp(:,m) = interp1(Ased(n).sigmacrop(1:index,m),Ased(n).ACvel(1:index,m),sigint);
        end
    end
    Ased(n).velinterp(Ased(n).velinterp == 0) = nan;
    Ased(n).StAvelprof = nan;
    Ased(n).StAvelprof = nanmean(Ased(n).velinterp,2);
    Ased(n).Adepthprof(:,1) = nan;
    Ased(n).Adepthprof(:,1) = sigint.*Stationdepth;
    Ased(n).Adepthprof = Stationdepth - Ased(n).Adepthprof;
end

for p = 1:length(Ased);
    Ased(p).lnAdepthcm = log((Ased(p).Adepthprof).*100);
    Ased(p).lnAdepthcm(Ased(p).lnAdepthcm<5.5) = nan; % to avoid extrapolated data
end

%% Now plot
for s = 1:length(Ased);
    figure;
    plot(abs(Ased(s).StAvelprof),Ased(s).lnAdepthcm,'ko');
    str1 = datestr(Ased(s).LTmean);
    title(str1,'FontSize',22);
    xlabel('Water Velocity (cm/s)','FontSize',14);
    ylabel('ln(Z) (cm)','FontSize',14);
    ylim([5 7]);
    str2 = ('A-A #');
    str3 = num2str(s);
    str4 = [str2, str3];
    set(gca,'Color','white');
    %saveas(gcf,str4,'jpg');
    %saveas(gcf,str4,'fig');
end

%% Sort velocity by depth to easier index points to fit line
for p =1%:length(Ased)
    [dval(:,p),dind(:,p)] = sort(Ased(p).lnAdepthcm, 'ascend');
    vel_sort1(:,p) = Ased(p).StAvelprof(dind(:,p));
end
depth_sort = dval;
vel_sort = vel_sort1;
% vel_sort1(vel_sort1 == inf) = nan;
% dval(isnan(vel_sort1)) = nan;
% vel_sort = relocateNaN(vel_sort1);
% depth_sort = relocateNaN(dval);
vel_sort = abs(vel_sort);

%% Fit Line to Bottom Boundary Layer
%First, look at graphs and choose what points to use
%U* = poly(:,1)
poly(1,:)=polyfit(vel_sort(1:7,1),depth_sort(1:7,1),1);
poly(2,:)=polyfit(vel_sort(1:11,2),depth_sort(1:11,2),1);
poly(3,:)=polyfit(vel_sort(1:8,3),depth_sort(1:8,3),1);
poly(4,:)=polyfit(vel_sort(1:8,4),depth_sort(1:8,4),1);
% poly(5,:)=polyfit(vel_sort(1:3,5),depth_sort(1:3,5),1);
% poly(6,:)=polyfit(vel_sort(1:6,6),depth_sort(1:6,6),1);
poly(7,:)=polyfit(vel_sort(2:5,7),depth_sort(2:5,7),1);
% poly(8,:)=polyfit(vel_sort(1:5,8),depth_sort(1:5,8),1);
% poly(9,:)=polyfit(vel_sort(1:8,9),depth_sort(1:8,9),1);
poly(10,:)=polyfit(vel_sort(1:13,10),depth_sort(1:13,10),1);
poly(11,:)=polyfit(vel_sort(1:4,11),depth_sort(1:4,11),1);
poly(12,:)=polyfit(vel_sort(1:28,12),depth_sort(1:28,12),1);
poly(13,:)=polyfit(vel_sort(1:28,13),depth_sort(1:28,13),1);
poly(14,:)=polyfit(vel_sort(1:10,14),depth_sort(1:10,14),1);
poly(15,:)=polyfit(vel_sort(1:9,15),depth_sort(1:9,15),1);
poly(16,:)=polyfit(vel_sort(1:13,16),depth_sort(1:13,16),1);
poly(17,:)=polyfit(vel_sort(1:3,17),depth_sort(1:3,17),1);
poly(18,:)=polyfit(vel_sort(1:3,18),depth_sort(1:3,18),1);
% poly(19,:)=polyfit(vel_sort(1:3,19),depth_sort(1:3,19),1);
poly(20,:)=polyfit(vel_sort(1:7,20),depth_sort(1:7,20),1);
poly(21,:)=polyfit(vel_sort(2:9,21),depth_sort(2:9,21),1);
poly(22,:)=polyfit(vel_sort(2:6,22),depth_sort(2:6,22),1);
poly(23,:)=polyfit(vel_sort(1:4,23),depth_sort(1:4,23),1);
poly(24,:)=polyfit(vel_sort(1:7,24),depth_sort(1:7,24),1);
uint=0:0.5:200;
polyv(1,:)=polyval(poly(1,:),uint);
polyv(2,:)=polyval(poly(2,:),uint);
polyv(3,:)=polyval(poly(3,:),uint);
polyv(4,:)=polyval(poly(4,:),uint);
polyv(5,:)=polyval(poly(5,:),uint);
polyv(6,:)=polyval(poly(6,:),uint);
polyv(7,:)=polyval(poly(7,:),uint);
polyv(8,:)=polyval(poly(8,:),uint);
polyv(9,:)=polyval(poly(9,:),uint);
polyv(10,:)=polyval(poly(10,:),uint);
polyv(11,:)=polyval(poly(11,:),uint);
polyv(12,:)=polyval(poly(12,:),uint);
polyv(13,:)=polyval(poly(13,:),uint);
polyv(14,:)=polyval(poly(14,:),uint);
polyv(15,:)=polyval(poly(15,:),uint);
polyv(16,:)=polyval(poly(16,:),uint);
polyv(17,:)=polyval(poly(17,:),uint);
polyv(18,:)=polyval(poly(18,:),uint);
polyv(19,:)=polyval(poly(19,:),uint);
polyv(20,:)=polyval(poly(20,:),uint);
polyv(21,:)=polyval(poly(21,:),uint);
polyv(22,:)=polyval(poly(22,:),uint);
polyv(23,:)=polyval(poly(23,:),uint);
polyv(24,:)=polyval(poly(24,:),uint);
%% Now plot with trendlines
for s = 1%[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24];
    plot(abs(Ased(s).StAvelprof),Ased(s).lnAdepthcm,'ko',uint,polyv(s,:),'r-');
    hold on
    if s == 1;
        plot(vel_sort(1:7,s),depth_sort(1:7,s),'ro'); 
    elseif s == 2;
        plot(vel_sort(1:11,s),depth_sort(1:11,s),'ro'); 
    elseif s == 3;
        plot(vel_sort(1:8,s),depth_sort(1:8,s),'ro');     
    elseif s == 4;
        plot(vel_sort(1:8,s),depth_sort(1:8,s),'ro');   
%     elseif s == 5;
%         plot(vel_sort(1:3,s),depth_sort(1:3,s),'ro');
%     elseif s == 6;
%         plot(vel_sort(1:6,s),depth_sort(1:6,s),'ro');
    elseif s == 7;
        plot(vel_sort(2:5,s),depth_sort(2:5,s),'ro');
%     elseif s == 8;
%         plot(vel_sort(1:5,s),depth_sort(1:5,s),'ro');  
%     elseif s == 9;
%         plot(vel_sort(1:8,s),depth_sort(1:8,s),'ro');       
    elseif s == 10;
        plot(vel_sort(1:13,s),depth_sort(1:13,s),'ro');
    elseif s == 11;
        plot(vel_sort(1:4,s),depth_sort(1:4,s),'ro');
    elseif s == 12;
        plot(vel_sort(1:28,s),depth_sort(1:28,s),'ro');    
    elseif s == 13;
        plot(vel_sort(1:28,s),depth_sort(1:28,s),'ro'); 
    elseif s == 14;
        plot(vel_sort(1:10,s),depth_sort(1:10,s),'ro');   
    elseif s == 15;
        plot(vel_sort(1:9,s),depth_sort(1:9,s),'ro');
    elseif s == 16;
        plot(vel_sort(1:13,s),depth_sort(1:13,s),'ro');
    elseif s == 17;
        plot(vel_sort(1:3,s),depth_sort(1:3,s),'ro');
    elseif s == 18;
        plot(vel_sort(1:3,s),depth_sort(1:3,s),'ro');  
%     elseif s == 19;
%         plot(vel_sort(1:7,s),depth_sort(1:7,s),'ro');       
    elseif s == 20;
        plot(vel_sort(1:7,s),depth_sort(1:7,s),'ro');
    elseif s == 21;
        plot(vel_sort(2:9,s),depth_sort(2:9,s),'ro');
    elseif s == 22;
        plot(vel_sort(2:6,s),depth_sort(2:6,s),'ro');    
    elseif s == 23;
        plot(vel_sort(1:4,s),depth_sort(1:4,s),'ro'); 
    elseif s == 24;
        plot(vel_sort(1:7,s),depth_sort(1:7,s),'ro');     
    end
    set(gca,'Color','white');
    str1 = datestr(Ased(s).LTmean);
    title(str1,'FontSize',22);
    xlabel('Water Velocity (cm/s)','FontSize',14);
    ylabel('ln(Z) (cm)','FontSize',14);
    ylim([3 7.5]);
    str2 = ('BBL A-A #');
    str3 = num2str(s);
    str4 = [str2, str3];
%     str5 = num2str(poly(s,1));
%     str6 = num2str(exp((poly(s,2))));
%     str7 = ('U* = ');
%     str8 = ('Zo = ');
%     str9 = [str7, str5, str8, str6];
    %saveas(gcf,str4,'jpg');
    %saveas(gcf,str4,'fig');
    hold off
end


%% Save U* and Zo
for n = 1%:length(Ased);
    Ased(n).UStarA = nan;
    Ased(n).UStarA = 0.41/poly(n,1)
    Ased(n).lnZoA = poly(n,2)
    Ased(n).ZoA = exp(poly(n,2))
    if Ased(n).ZoA ==1;
        Ased(n).ZoA = nan;
    end
end

%% Calculate Umean for several reference heights and find Cd
for n = 1:length(Ased);
[val_depth350,i_depth350] = min(abs((Ased(1).Adepthprof)-3.5));
Ased(n).UmeanA350 = abs(Ased(n).StAvelprof(i_depth350));        
        
[val_depth300,i_depth300] = min(abs((Ased(1).Adepthprof)-3));
Ased(n).UmeanA300 = abs(Ased(n).StAvelprof(i_depth300));        
        
[val_depth250,i_depth250] = min(abs((Ased(1).Adepthprof)-2.5));
Ased(n).UmeanA250 = abs(Ased(n).StAvelprof(i_depth250));

[val_depth200,i_depth200] = min(abs((Ased(1).Adepthprof)-2));
Ased(n).UmeanA200 = abs(Ased(n).StAvelprof(i_depth200));
end

for s = 1:length(Ased);
    UstarA(s,1) = Ased(s).UStarA;
    UmeanA350(s,1) = Ased(s).UmeanA350;
    UmeanA300(s,1) = Ased(s).UmeanA300;
    UmeanA250(s,1) = Ased(s).UmeanA250;
    UmeanA200(s,1) = Ased(s).UmeanA200;
end

UstarA(UstarA == inf) = nan;
polyindex = find(~isnan(UstarA));

% Cd is slope from U*^2 against Umean^2
CdPoly350 = polyfit((UmeanA350(polyindex)).^2,UstarA(polyindex).^2,1);
CdPoly300 = polyfit((UmeanA300(polyindex)).^2,UstarA(polyindex).^2,1);
CdPoly250 = polyfit((UmeanA250(polyindex)).^2,UstarA(polyindex).^2,1);
CdPoly200 = polyfit((UmeanA200(polyindex)).^2,UstarA(polyindex).^2,1);

Cd350 = CdPoly350(1);
Cd300 = CdPoly300(1);
Cd250 = CdPoly250(1);
Cd200 = CdPoly200(1);

for n = 1:length(Ased);
    Ased(n).Cd350A = CdPoly350(1);
    Ased(n).Cd300A = CdPoly300(1);
    Ased(n).Cd250A = CdPoly250(1);
    Ased(n).Cd200A = CdPoly200(1);
end

plot((UmeanA350(polyindex)).^2,UstarA(polyindex).^2, 'ro', (UmeanA300(polyindex)).^2, UstarA(polyindex).^2,'bo',(UmeanA250(polyindex)).^2, UstarA(polyindex).^2,'ko');
hold on
title ('A-A2014 U* vs Umean', 'FontSize', 20);
xlabel ('Umean^2 (cm/s)^2', 'FontSize', 16);
ylabel ('U*^2 (cm/s)^2', 'FontSize', 16);
str1 = ['Cd350 = ',num2str(Cd350)];
str2 = ['Cd300 = ',num2str(Cd300)];
str3 = ['Cd250 = ',num2str(Cd250)];
legend(str1, str2, str3, 'Location', 'best');
saveas(gcf,'A-ACd2014','jpg');
saveas(gcf,'A-ACd2014','fig');

%% Use Cd calculated from Law of the Wall to find U* for all passes
% Tau = density * Cd250 * Umean250^2 = density * Ustar^2
% Cd250 * Umean250^2 = Ustar^2
% Ustar = sqrt (Cd250 * Umean250^2)

for n = 1:length(Ased);
    if n ~=3;
    Ased(n).UstarAQSL = sqrt(Ased(n).Cd250A .* (Ased(n).UmeanA250).^2);
        if Ased(n).UStarA == inf;
            Ased(n).UStarA = NaN;
        end
    end
end

for n = 1:length(Ased);
    if n~=3;
    UstarAQSL(n,1) = Ased(n).UstarAQSL;
    UstarA(n,1) = Ased(n).UStarA;
    time(n,1) = Ased(n).LTmean;
    end
end

plot(time,UstarA, 'ro');
hold on
plot(time, UstarAQSL, 'ko');
title('Transect-A(2014) U*', 'FontSize', 20);
datetick('x', 'HH:MM');
xlabel('Local Time', 'FontSize', 16);
ylabel('U* (cm/s)', 'FontSize', 16);
legend('Law of the Wall', 'Quadratic Stress Law');
hFig = figure(1);
set(hFig, 'Position', [500 500 1000 500])
saveas(gcf,'A-AUstar2014','jpg');
saveas(gcf,'A-AUstar2014','fig');






