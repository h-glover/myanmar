clear all,close all,clc

% compare elevation to water level
load('C:\GLOVER\data\myanmar\SurveyPts2019.mat')
load('C:\GLOVER\output\myanmar\longterminst\BogaleRiverInstruments.mat')
br = BogaleRiver; clear BogaleRiver

site_name={'croc station';'cyclone';'freda';'monkey';'cyclone';'rice'};

msl_hf = nanmean(br.MeinmahlaDepth(132900:144100));
br.MeinmahlaWaterLevel=br.MeinmahlaDepth-msl_hf; 
hhw_hf = max(br.MeinmahlaWaterLevel(132900:144100));
llw_hf = min(br.MeinmahlaWaterLevel(132900:144100));
idx = 133387:143954;
coefs = ut_solv(br.datenum(idx)',br.FredaDepth(idx),[],16,'auto');
pred = ut_reconstr(br.datenum(idx),coefs);
msl_hf_freda = nanmean(pred);
%%
% % calculate the mean and the mean higher high water levels for the HF
% % record:
% mhhw_hf = br.MeinmahlaDepth(132900:132899+(72*24*60/10));
% mhhw_hf = reshape(mhhw_hf,12*60/10,[]);
% mhhw_hf = max(mhhw_hf,[],2);
% mhhw_hf = nanmean(mhhw_hf);

% pull out all data for water levels
flds = fieldnames(survey);
for jj=1:length(flds)
    wl.(flds{jj})=survey.(flds{jj})([20,58,75,101,119,147,185:187]);
end
wltimes = {'9/30/2019 11:40 AM';'9/30/2019 12:00 PM';'9/30/2019 8:20 AM';...
    '9/29/2019 1:30 PM';'9/28/2019 5:40 PM';'9/28/2019 10:40 AM';...
    '9/28/2019 8:20 AM';'9/28/2019 8:20 AM';'9/28/2019 8:30 AM'};
wl.datenum = datenum(wltimes);
[wl.datenum,srt_idx]=sort(wl.datenum,'ascend');
for jj=1:length(flds)
    wl.(flds{jj})=wl.(flds{jj})(srt_idx);
end
% dz = water level during survey + survey elev to that datum
water_edge_elev = abs(wl.Zcorr(4));
water_edge2msl = br.MeinmahlaWaterLevel(br.datenum==wl.datenum(4));
dz = water_edge_elev + water_edge2msl;
% correct all survey elevations to the HF mean water level
survey.Zcorr_wl = survey.Zcorr + dz;
wl.Zcorr_wl = wl.Zcorr + dz;

survey.ID(survey.ID==5)=2;
for jj=[1:4,6]
    figure;
    if jj==1 || jj==2 || jj==5
        xx = survey.Lat(survey.ID==jj & survey.type==0);
    elseif jj==3 || jj==4 || jj==6
        xx = survey.Lon(survey.ID==jj & survey.type==0);
    end
    plot(xx,survey.Zcorr_wl(survey.ID==jj & survey.type==0),'*'),hold on
    R=refline(0,0);R.Color='k';
    R=refline(0,hhw_hf);R.Color='r';R=refline(0,llw_hf);R.Color='r';
    ylim([-4 4])
    title(site_name{jj})
end


% figure;plot(survey.Lat,survey.Zcorr_wl,'*')
% xlabel('latitude'),ylabel('elevation wrt msl at hf')
% ylim([-4 4]),refline(0,0)

% only plot actual land surface values - not structures
figure;
idx = find(survey.type==0);
scatter(survey.Lon(idx),survey.Lat(idx),[],survey.Zcorr_wl(idx),'filled')
hold on
idx = find(survey.type==2);
scatter(survey.Lon(idx),survey.Lat(idx),[],survey.Zcorr_wl(idx),'*')
colorbar,caxis([-3 3]),colormap(cmocean('balance'))

figure;
plot(survey.ID,survey.Z,'k*'),hold on
plot(survey.ID,survey.Zcorr,'r*')

%% figure(1),refline(0,survey.Z(151)),refline(0,survey.Z(152))
figure(3),refline(0,survey.Z(102))

figure;idx = find(survey.ID==3);
scatter([0:length(idx)-1],survey.Z(idx),[],survey.Lon(idx),'*')
colorbar%,caxis([-3 3])
% hold on,scatter(survey.Lon(102),survey.Lat(102),'ko')
refline(0,0)

figure;
for jj=1:6
    idx = find(survey.ID==jj);
    plot(survey.Lat(idx),survey.Z(idx),'*'),hold on
end
legend
R=refline(0,0);R.Color='k';
% R=refline(0,hhw_hf);R.Color='r';R=refline(0,llw_hf);R.Color='r';
% refline(0,survey.Z(121)),refline(0,survey.Z(59))
refline(0,survey.Z(102))

% plot of data overlaps
load('C:\GLOVER\output\myanmar\aqd\aqd_sep19_hc.mat'),hc=aqd;
load('C:\GLOVER\output\myanmar\aqd\aqd_sep19_lc.mat'),lc=aqd;
load('C:\GLOVER\output\myanmar\aqd\aqd_sep19_ag.mat'),ag=aqd;

figure;
% subplot(121)
plot(hc.time,hc.depth,'k--'),hold on
plot(lc.time,lc.depth,'k-.')
plot(ag.time,ag.depth,'k:')
plot(br.datenum,br.FredaDepth,'b--')
plot(br.datenum,br.MeinmahlaWaterLevel,'-b')
scatter(wl.datenum,wl.Z,[],wl.Lat,'*')
xlim([datenum('9/27/2019') datenum('10/1/2019')])
datetick('x','dd','keeplimits')
legend({'aqd hc','aqd lc','aqd ag','Freda dock','Cyclone Dock','wl survey pts'})
colorbar
% subplot(122),scatter(wl.Lon,wl.Lat,[],[1:8]','*'),colorbar,caxis([1 7])

% ref croc hc to wl 1&2
offset = hc.depth(hc.time==wl.datenum(1)) - mean([wl.Z(1),wl.Z(2)]);
hc.waterelev = hc.depth - offset;
plot(hc.time,hc.waterelev,'r--')

% ref cyclone lc to wl 4
offset = lc.depth(lc.time==wl.datenum(4)) - wl.Z(4);
lc.waterelev = lc.depth - offset;
plot(lc.time,lc.waterelev,'r-.')

% ref cyclone ag to wl 8
t = wl.datenum(8) - datenum(0,0,0,12,25,0);
offset = ag.depth(ag.time==t) - wl.Z(8);
ag.waterelev = ag.depth - offset;
plot(ag.time,ag.waterelev,'r:')

% ref freda to wl 5


figure;
plot(br.datenum,br.MeinmahlaWaterLevel,'k'),hold on
plot(br.datenum,br.FredaDepth,'b-')
plot(br.datenum,br.ut_FredaDepth,'b:')
refline(0,0)



%% survey map data
clear all,close all,clc
cd 'C:\GLOVER\data\myanmar\'
fid=fopen('SurveyPts2019.csv');
T=textscan(fid,'%*f %*f %*f %*f %s %f %f %f %f %f %f %f %f',...
    'Delimiter',',','HeaderLines',1);
fclose(fid);


survey.site=T{1};
survey.ID=T{8};
survey.Lat=T{3};
survey.Lon=T{2};
survey.Z=T{4};
survey.type=T{9};

% compare ref to sanity checks and use that offset to fix???
survey.Zcorr = survey.Z;
% Site1: Ref=151; San=152;
survey.Zcorr(152:end) = survey.Z(152:end)+(survey.Z(151)-survey.Z(152));
survey.Zcorr([182,183,185,186])=2+survey.Zcorr([182,183,185,186]);
% Site2: Ref=121; San=[122:123];San_basetower=[124:125];
survey.Zcorr(122:150) = survey.Z(122:150)+(survey.Z(121)-survey.Z(122));
% Site3: Ref=102; San=103;
% survey.Zcorr(102:120) = survey.Z(102:120)+(survey.Z(102)-survey.Z(103));
% Site4: Ref=81; San=82;
survey.Zcorr(82:93) = survey.Z(82:93)+(survey.Z(81)-survey.Z(82));
% Site4: Ref=94; San=99;
survey.Zcorr(95:101) = survey.Z(95:101)+(survey.Z(94)-survey.Z(99));
survey.Zcorr(94:101) = survey.Zcorr(94:101)+(survey.Z(93)-survey.Z(94));
% Site5: Ref=59; San=76;
% survey.Zcorr(60:79) = survey.Z(60:79)+(survey.Z(59)-survey.Z(76));
% Site6: Ref=1; San=2;
survey.Zcorr(2:58) = survey.Z(2:58)+(survey.Z(1)-survey.Z(2));
% type of measurement: 0=land, 1=structure, 2=maybe?


figure;
plot(survey.ID,survey.Z,'k*'),hold on
plot(survey.ID,survey.Zcorr,'r*')


% plot only land surface values?

save('SurveyPts2019','survey')