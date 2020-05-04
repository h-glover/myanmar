clear all,close all,clc

ff=dir('*_FluxDecomp2.mat');
for jj=1:length(ff)
    load(ff(jj).name)
    di(1,jj)=round(fluxdecomp.Di0,-1);
    dierr(1,jj)=fluxdecomp.rmse.Di0;
    F_comp(jj)=fluxdecomp.Fs0;
    Fe_comp(jj)=fluxdecomp.Fe;
    Fr_comp(jj)=fluxdecomp.Fr;
    Ft_comp(jj)=fluxdecomp.Ft;
    Ferr_comp(jj)=fluxdecomp.rmse.Fs;
    Fdiff_comp(jj)=F_comp(jj)-(Fe_comp(jj)+Fr_comp(jj));
end
lbl_BR={'LF, Neap','LF, Spr','HF, Neap','HF, Spr'};
lbl_YR={'LF, Neap','HF, Neap'};
M=[F_comp;Fr_comp;Fe_comp;Ft_comp;Ferr_comp];

figure;
subplot(211)
bar(M(:,1:4)'),hold on
ax=gca; ax.XTickLabels=lbl_BR;
ylabel('Sediment Flux (t/s)')
title('Bogale River')
ylim([-0.5 0.75])

subplot(212)
bar(M(:,5:6)')
ax=gca; ax.XTickLabels=lbl_YR;
ylabel('Sediment Flux (t/s)')%,legend('F_{TOTAL}','F_R','F_E','F_T')
title('Yangon River')
ylim([-2.5 3.5])


%%
clear all,close all,clc

ff=dir('*_SaltFluxDecomp3.mat');
for jj=1:length(ff)
    load(ff(jj).name)
    F_comp(jj)=fluxdecomp.Fs0;
    Fe_comp(jj)=fluxdecomp.Fe;
    Fr_comp(jj)=fluxdecomp.Fr;
    Ft_comp(jj)=fluxdecomp.Ft;
    Ferr_comp(jj)=fluxdecomp.rmse.Fs;
    Fdiff_comp(jj)=F_comp(jj)-(Fe_comp(jj)+Fr_comp(jj));
end
lbl_BR={'LF, Neap','LF, Spr'};
lbl_YR={'LF, Neap','LF, Spring'};
M=[F_comp;Fr_comp;Fe_comp;Ft_comp;Ferr_comp;Fdiff_comp];
M(:,4)=NaN;
%
figure;
subplot(211)
bar(M(:,1:2)'),hold on
ax=gca; ax.XTickLabels=lbl_BR;
title('Bogale River')
% ylim([-6.5 5])

subplot(212)
bar(M(:,3:4)')
ax=gca; ax.XTickLabels=lbl_YR;
title('Yangon River')
% ylim([-2.5 5])


%%

clear all,close all,clc

ff=dir('*_FluxDecomp3.mat');
for jj=1:length(ff)
    load(ff(jj).name)
    F_comp(jj)=fluxdecomp.Fs0;
    Fe_comp(jj)=fluxdecomp.Fe;
    Fr_comp(jj)=fluxdecomp.Fr;
    Ft_comp(jj)=fluxdecomp.Ft;
end
lbl_BR={'LF, Neap','LF, Spr','HF, Neap','HF, Spr'};
lbl_YR={'LF, Neap','HF, Neap'};
M=[F_comp;Fr_comp;Fe_comp;Ft_comp];
%
figure;
subplot(221)
bar(M(:,1:4)'),hold on
ax=gca; ax.XTickLabels=lbl_BR;
ylabel('Sediment Flux (t/s)')
title('Bogale River')
ylim([-0.5 0.75])

subplot(223)
bar(M(:,5:6)')
ax=gca; ax.XTickLabels=lbl_YR;
ylabel('Sediment Flux (t/s)')%,legend('F_{TOTAL}','F_R','F_E','F_T')
title('Yangon River')
ylim([-2.5 3.5])
%%
clear all
ff=dir('*SaltFluxDecomp2.mat');
for jj=1:length(ff)
    load(ff(jj).name)
    F_comp(jj)=fluxdecomp.Fs0;
    Fe_comp(jj)=fluxdecomp.Fe;
    Fr_comp(jj)=fluxdecomp.Fr;
    Ft_comp(jj)=F_comp(jj)-(Fe_comp(jj)+Fr_comp(jj));
    %Ft_comp(jj)=fluxdecomp.Ft;
    disch(jj)=fluxdecomp.U0*fluxdecomp.A0;
end

lbl_BR={'LF, Neap','LF, Spr'};
lbl_YR={'LF, Neap','LF, Spr'};
M=[F_comp;Fr_comp;Fe_comp;Ft_comp];


subplot(222)
bar(M(:,1:2)')
ax=gca; ax.XTickLabels=lbl_BR;
ylabel('Salt Flux (t/s)')
title('Bogale River')
ylim([-6.5 4])

subplot(224)
bar([M(:,3)';[NaN NaN NaN NaN]])
ax=gca; ax.XTickLabels=lbl_YR;
ylabel('Salt Flux (t/s)'),legend('F_{TOTAL}','F_R','F_E','F_T')
title('Yangon River')
ylim([-6.5 4])

%%
load('GrainSizeData_Lab4.mat')
close all
figure;
ttl={'Stn 5','Stn 23','Stn24'};
for jj=1:3
    
    subplot(1,3,jj)
    semilogx(x(:,1),x(:,jj+1))
    xlim([0.0001 1]),ylim([0 70]),title(ttl{jj}) 
end
