
clear all,close all,clc

cd C:\GLOVER\output\myanmar\pb210
% plot all except: 323(1), 324(2), 347(7), 432(10); 524(16) + 463(12) = 2b; 
F=dir('*ex.mat');
d_bot = [100,100,30,30,40,35,30,60,50,60,50,50,60,100,30,60,25,20,40,100];
idx=[3:13,17:19];
    
figure; cc=1;
for jj=idx
    
    % load the data and compile all the core intervals into 1 vector for
    % each variable
    load(F(jj).name)
    depth = vertcat(pb210.midinterval);
    act = vertcat(pb210.actatcoll_saltmudcor);
    act(act<0)=NaN;
    act_err = vertcat(pb210.tot_err)./2;
    depth_err = vertcat(pb210.intervalsize)./2;
    
    % remove a few select bad surface values
    if jj==3 || jj==8 || jj==10 || jj==12
        act(3)=NaN;
    elseif jj==3% || jj==14 || jj==15
        act(5)=NaN;
    elseif jj==9 || jj==11 || jj==15 || jj==16
        act(1)=NaN;
    end

% plot the raw data
%     figure(fnum(jj));
%     subplot(1,5,snum(jj))
    subplot(2,7,cc)
    errorbar(act,depth,depth_err,depth_err,act_err,act_err,'.')
    axis ij,hold on,
    ax=gca; ax.XScale='log';
    xlim([0.1 5])
    ylim([0 100])
    
    % Code to define the data range over which to fit the line
    xx=act(depth<d_bot(jj) & ~isnan(act));
    yy=depth(depth<d_bot(jj) &~isnan(act));
    % Linear Regression
    b=polyfit(log10(xx),yy,1);
    x_fit=[min(xx),max(xx)];
    y_fit = polyval(b,log10(x_fit));
    
    % plot the fit line
    % plot(xx,yy,'ro')
    semilogx(x_fit,y_fit,'k-')
    
    % Text
    [~,maxy]=max(y_fit);
    [~,miny]=min(y_fit);
    s =(0.0311*range(y_fit))/(log(x_fit(2)/x_fit(1)));
    s = round(s,2);
    % calculate the r2 fit for the accumulation rate line
    y_fit = polyval(b,log10(xx));
    r = corrcoef(yy,y_fit);
    R2(jj) = round(r(2,1)^2,2);
    
    title([F(jj).name(1:5),'\newlineA=',num2str(s),' cm/y, ',num2str(R2(jj))])
    
    cc=cc+1;
end


%% grain size for all plots above:
clear all,close all,clc
load('C:\GLOVER\output\myanmar\grainsize\BR_GrainSize_Cores2.mat')
cores.br346 = cores.br346([1,3:9]);
cores.br427 = cores.br427([1,2,4,6:8]);
cores.br520 = cores.br520([2:8]);
cores.br521 = cores.br521([1,2,4:10]);
cores.br522 = cores.br522([2:9]);
corename= {'br323','br324','br326','br327','br335','br346','br360',...
    'br427','br428','br429','br430','br431','br346','br462',...
    'br520','br521','br522'};
% 347(3), 432(10) 524(16) + 463(12) = 2b; 

figure;
for jj=1:14
% plot grain size
frac = cores.(corename{jj});
frac = vertcat(frac.frac)/100;

subplot(2,7,jj)
barh(depth,frac,'stacked')
axis ij
ylim([0 100]), xlim([0 1])

end
legend({'clay','silt','sand'})

