function [S]= conduc2sali(C,T,P)
% _________________________________________________
% This script computes the water salinity from pressure and
% conductivity of water data from Fofonoff and Millard (1983)
% equation. User can compute a single value or an array 
% of values. 
% 
% The array with values must have the next order: Column 1,
% water conductivity; Column 2, water temperature; and Column 3
% water pressure.
%
% Input:
%          C = Water conductivity, mS/cm 
%          T = water temperature, C
%          P = water pressure, db.
%
% Output:
%          salin1 = water salinity, PSS-78 units
%
% Examples:
%          sal = conduc2sali(30,10,1000)
%          data = [ 30 20 0; 40 20 1000; 50 20 2000];
%          sal = conduc2sali(data)
%
% Gabriel Ruiz Mtz.
% Version 2
% 09/10/10
% Rev. 2017.08.08: Fix the conversion from S/m to mS/cm
% _________________________________________________

mi_in = 1;
ma_in = 3;
narginchk(mi_in,ma_in)

%% If we have an array...
if nargin == 1
    du = C;
    C = du(:,1);
    T = du(:,2);
    P = du(:,3);
    clear mi_in ma_in du;
elseif nargin == 2
    error('Only 2 inputs arguments were used, when the function requires 3');
end

%% Computing the salinity
S = salinity(C,T,P);

%% Ploting the single value
% You can turn off the line 50 to 65 if you do not want the plot.
% if nargin == 3
%    gp= 0:10:10000;
%    gc = ones(length(gp),1);
%    gt = ones(length(gp),1);
%    gc = gc.*C;
%    gt = gt.*T;
%    gs = salinity(gc,gt,gp');
% %    plot(gp,gs,'LineWidth',2);
% %    xlabel('Water pressure, [dB]');
% %    ylabel('Water salinity, [PSS-78]');
% %    hold on;
% %    h1 = plot(P,S,'r*','MarkerSize',8);
% %    legend(h1,{horzcat('P= ',num2str(P),' dB, T = ',num2str(T),' °C, C = ',num2str(C),' mS/cm')});
% %    set(gca,'Xgrid','on','Ygrid','on');
% end

return

function [S] = salinity(C,T,P)
% ________________________________________________________
% Subfunction to compute the water salinity from values
% of water conductivity, temperature and pressure.
%
% References:
% Fofonoff, P. and Millard, R.C. (1983). Algorithms for 
%    computation of fundamental properties of seawater,  
%    Unesco Technical Papers in Marine Sci., 44, 58 pp.
% UNESCO. (1981). Background papers and supporting data on
%    the practical salinity, 1978. Unesco Technical Papers 
%    in Marine Sci., 37, 144 pp.
% Wagner, R.J., Boulger, R.W., Jr., Oblinger, C.J., y Smith, 
%    B.A., 2006, Guidelines and standard procedures for continuous
%    water-quality monitorsStation operation, record computation, 
%    and data reporting: U.S. Geological Survey Techniques
%    and Methods 1D3, 51 p.; access 2006-04-10 
%    en http://pubs.water.usgs.gov/tm1d3
% http://www.salinometry.com/welcome/
%
% Gabriel Ruiz Mtz.
% Version 2
% ________________________________________________________

%% Computing the conductivity ratio
% $ R = /frac{C(S,T_{68},P)}{C(35,15_{68},0)} $
% where $ C(35,15_{68},0) $ = 42.914 mS/cm = 4.2914 S/m
cnd = C/42.914;
R = cnd;
     
%% Computing rt(t)
if (min(T) >= -2) && (max(T) <= 35)
    c0 =  0.6766097;
	c1 =  2.00564e-2;
	c2 =  1.104259e-4;
	c3 =  -6.9698e-7;
	c4 =  1.0031e-9;
    RT35 = ( ( (c3+c4.*T).*T+c2).*T+c1).*T+c0;
else
	error('Temperature out of range');
end

%% Computing Rp(S,t,p)
d1 = 3.426e-2;
d2 = 4.464e-4;
d3 = 4.215e-1;
d4 = -3.107e-3;
e1 = 2.070e-5;
e2 = -6.370e-10;
e3 = 3.989e-15;
RP = 1+(P.*(e1+e2.*P+(e3.*P.^2)))./(1+(d1.*T)+(d2.*T.^2)+(d3+(d4.*T)).*R);
     
%% Computing Rt(S,t)
RT = R./(RP.*RT35);
	 
%% Computing R
XR =sqrt(RT);
XT = T - 15;
a0 = 0.0080;
a1 = -0.1692;
a2 = 25.3851;
a3 = 14.0941;
a4 = -7.0261;
a5 = 2.7081;
b0 =  0.0005;
b1 = -0.0056;
b2 = -0.0066;
b3 = -0.0375;
b4 =  0.0636;
b5 = -0.0144;
k  =  0.0162;
DSAL = (XT./(1+k.*XT)).*(b0+(b1+(b2+(b3+(b4+(b5.*XR)).*XR).*XR).*XR).*XR);
SAL = (((((a5.*XR)+a4).*XR+a3).*XR+a2).*XR+a1).*XR+a0;
S = SAL + DSAL;
return
