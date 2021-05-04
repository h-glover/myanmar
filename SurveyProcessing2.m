clear all,close all,clc

% compare elevation to water level
load('C:\GLOVER\data\myanmar\SurveyPts2019.mat')
load('C:\GLOVER\output\myanmar\longterminst\BogaleRiverInstruments.mat')
br = BogaleRiver; clear BogaleRiver

site_name={'croc station';'cyclone';'freda';'monkey';'cyclone';'rice'};

msl_hf_m = nanmean(br.MeinmahlaDepth(132900:144100));
hhw_hf = max(br.MeinmahlaWaterLevel(132900:144100));
llw_hf = min(br.MeinmahlaWaterLevel(132900:144100));

% pull out all data for water levels
flds = fieldnames(survey);
for jj=1:length(flds)
    wl.(flds{jj})=survey.(flds{jj})([20,58,75,101,119,147,185:187]);
end
wltimes = {'9/30/2019 11:40 AM';'9/30/2019 12:00 PM';'9/30/2019 8:20 AM';...
    '9/29/2019 1:30 PM';'9/28/2019 5:40 PM';'9/28/2019 10:40 AM';...
    '9/28/2019 8:20 AM';'9/28/2019 8:20 AM';'9/27/2019 7:40 AM'};
wl.datenum = datenum(wltimes);
[wl.datenum,srt_idx]=sort(wl.datenum,'ascend');
for jj=1:length(flds)
    wl.(flds{jj})=wl.(flds{jj})(srt_idx);
end

figure;
yyaxis right,plot(wl.datenum,wl.Zcorr,'o'),hold on
plot(wl.datenum(4),wl.Zcorr(4),'*')
plot(wl.datenum(1),wl.Zcorr(1),'*')
yyaxis left
plot(br.datenum,br.MeinmahlaWaterLevel,'k'),hold on
plot(br.datenum,br.FredaWaterLevel,'r')


% high tide is 30 minutes later at freda than at cyclone shelter
% so, to match the wl at croc station you would need to subtrct 30 min from
% timing of water level measurements in the survey
%Freda dz = water level during survey + survey elev to that datum
water_edge_elev = abs(wl.Zcorr(1));
water_edge2msl = br.FredaWaterLevel(br.datenum==(wl.datenum(1)));
dz_f = water_edge_elev + water_edge2msl;
% correct all survey elevations to the HF mean water level
survey.Zcorr_wl_f = survey.Zcorr + dz_f;
wl.Zcorr_wl_f = wl.Zcorr + dz_f;

%MMI dz = water level during survey + survey elev to that datum
water_edge_elev = abs(wl.Zcorr(4));
water_edge2msl = br.MeinmahlaWaterLevel(br.datenum==wl.datenum(4));
dz_m = water_edge_elev + water_edge2msl;
% correct all survey elevations to the HF mean water level
survey.Zcorr_wl_m = survey.Zcorr + dz_m;
wl.Zcorr_wl_m = wl.Zcorr + dz_m;

survey.Zcorr_wl = survey.Zcorr_wl_m;
survey.Zcorr_wl(survey.ID==1)=survey.Zcorr_wl_f(survey.ID==1);

survey.ID(survey.ID==5)=2;
for jj=[1:2,6]
    figure;
    if jj==1 || jj==2 || jj==5
        xx = survey.Lat(survey.ID==jj & survey.type==0);
    elseif jj==3 || jj==4 || jj==6
        xx = survey.Lon(survey.ID==jj & survey.type==0);
    end
    plot(xx,survey.Zcorr_wl(survey.ID==jj & survey.type==0),'k*'),hold on

    R=refline(0,0);R.Color='k';
    R=refline(0,hhw_hf);R.Color='r';R=refline(0,llw_hf);R.Color='r';
    ylim([-4 4])
    title(site_name{jj})
end

survey.msl_hf = msl_hf_m;
survey.hhw_hf = hhw_hf;
survey.llw_hf = llw_hf;

save('C:\GLOVER\output\myanmar\SurveyPts2019.mat','survey')

%% reference all the aquadopps to this 
clear all,close all,clc
cd C:\GLOVER\output\myanmar\
% compare elevation to water level
load('SurveyPts2019.mat')
load('longterminst\BogaleRiverInstruments.mat')
br = BogaleRiver; clear BogaleRiver

% sept2019 aquadopp data
load('aqd\aqd_sep19_hc.mat'),hc=aqd;
load('aqd\aqd_sep19_lc.mat'),lc=aqd;
load('aqd\aqd_sep19_ag.mat'),ag=aqd;

% pull out all data for water levels
flds = fieldnames(survey);
for jj=1:10
    wl.(flds{jj})=survey.(flds{jj})([20,58,75,101,119,147,185:187]);
end
wltimes = {'9/30/2019 11:40 AM';'9/30/2019 12:00 PM';'9/30/2019 8:20 AM';...
    '9/29/2019 1:30 PM';'9/28/2019 5:40 PM';'9/28/2019 10:40 AM';...
    '9/28/2019 8:20 AM';'9/28/2019 8:20 AM';'9/28/2019 8:10 AM'};
wl.datenum = datenum(wltimes);
[wl.datenum,srt_idx]=sort(wl.datenum,'ascend');
for jj=1:10
    wl.(flds{jj})=wl.(flds{jj})(srt_idx);
end

figure;
plot(hc.time,hc.depth,'k--'),hold on
plot(lc.time,lc.depth,'k-.')
plot(ag.time,ag.depth,'k:')
plot(br.datenum,br.FredaWaterLevel,'b--')
plot(br.datenum,br.MeinmahlaWaterLevel,'-b')
plot(wl.datenum,wl.Zcorr_wl,'o')
xlim([datenum('9/27/2019') datenum('10/1/2019')])
datetick('x','dd','keeplimits')
legend({'aqd hc','aqd lc','aqd ag','Freda dock','Cyclone Dock','wl survey pts'})

% pt1 is for the hc aqd:
hc.offset = wl.Zcorr_wl(1) - hc.depth(hc.time==wl.datenum(1));
hc.depth_elev=hc.depth+hc.offset;
hc.flat = nanmean(survey.Zcorr_wl(161:176));
plot(hc.time,hc.depth_elev,'g--')
% aqd=hc; save('aqd\aqd_sep19_hc.mat','aqd')
% 
ag.offset = wl.Zcorr_wl(8) - (ag.depth(end)+0.2);
ag.depth_elev=ag.depth+ag.offset;
ag.flat = nanmean(survey.Zcorr_wl([3:15,34:54]));
plot(ag.time,ag.depth_elev,'g:')
aqd=ag; save('aqd\aqd_sep19_ag.mat','aqd')

lc.offset = wl.Zcorr_wl(4) - lc.depth(1);
lc.depth_elev=lc.depth+lc.offset;
lc.flat= nanmean(survey.Zcorr_wl(126:139));
plot(lc.time,lc.depth_elev,'g-.')
% aqd=lc; save('aqd\aqd_sep19_lc.mat','aqd')
