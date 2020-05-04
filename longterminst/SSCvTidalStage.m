%% Comparing SSC to Tidal Stage for the Cyclone shelter channel

clear all,close all,clc
load('BogaleRiverInstruments.mat')

br=BogaleRiver; clear BogaleRiver

% % smooth SSC/SAL
% cidx=repmat(1/3,1,3);
% br.MeinmahlaSSC=conv(br.MeinmahlaSSC,cidx,'same');
% br.MeinmahlaSal=conv(br.MeinmahlaSal,cidx,'same');


% pull out 1 val per hour
flds=fieldnames(br);
tvec=br.datenum(1):datenum(0,0,0,1,0,0):br.datenum(end);
[~,inew,~]=intersect(br.datenum,tvec);
for jj=[1,3:length(flds)]
    br.(flds{jj})=br.(flds{jj})(inew);
end

% calculate water slope in m/hr: 
% wl(2) - wl (1) = m/hr
stg=(br.MeinmahlaDepth(2:end)-br.MeinmahlaDepth(1:end-1));
stg(end+1)=NaN;
stg(stg<-0.6 | stg>0.85)=NaN;
stg_ut=(br.ut_MeinmahlaDepth(2:end)-br.ut_MeinmahlaDepth(1:end-1));
stg_ut(end+1)=NaN;

% pull out HF and LF months:
dvec=datevec(br.datenum);
hf=find(dvec(:,2)>=7 & dvec(:,2)<11);
lf=find(dvec(:,2)<7 | dvec(:,2)>=11);

figure;
subplot(221)
scatter(stg(hf),br.MeinmahlaSSC(hf)),title('HF SSC')
subplot(222)
scatter(stg(lf),br.MeinmahlaSSC(lf)),title('LF SSC')
subplot(223)
scatter(stg(hf),br.MeinmahlaSal(hf)),title('HF Sal')
subplot(224)
scatter(stg(lf),br.MeinmahlaSal(lf)),title('LF Sal')


% reds=[ones(100,1),linspace(0,0.95,100)',linspace(0,0.95,100)'];
% blues=flipud([linspace(0,0.95,100)',linspace(0,0.95,100)',ones(100,1)]);
% C=[reds;blues];colormap(C)

reds=[0.9*ones(100,1),linspace(0.9,0,100)',linspace(0.9,0,100)'];
figure;
subplot(121)
scatter(stg_ut(hf),br.MeinmahlaSSC(hf),'k.'),hold on
scatter(stg_ut(hf),br.MeinmahlaSSC(hf),[],br.MeinmahlaSal(hf)),title('HF SSC')
colorbar,caxis([0 3])
subplot(122)
scatter(stg_ut(lf),br.MeinmahlaSSC(lf),'k.'),hold on
scatter(stg_ut(lf),br.MeinmahlaSSC(lf),[],br.MeinmahlaSal(lf)),title('LF SSC')
colorbar,caxis([1 15])
colormap(reds)