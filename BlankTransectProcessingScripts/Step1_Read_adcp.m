%% Read in ADCP data (Legacy ASCII text files from WinRiver).  Requires that you have the function adcp_ascii (D. Nowacki)

% Folder containing only text files of each transect
FileList=dir('*ASC.TXT');

for n = 1:length(FileList)
    disp(['processing ' num2str(n)]);
    
    adcp(n) = adcp_ascii(FileList(n).name);

end

% if you're satisfied with the matric produced, make sure to save it in
% your workspace!
save('test','adcp')
%%
close all

for jj=[1:3,20:23]%length(adcp)
    figure;
    L=1:length(adcp(jj).intensscale);
    pcolor(L,adcp(jj).z(:,1),adcp(jj).dir),shading flat
    colorbar,axis ij,caxis([1 360])

    %get the matlab time for each transect and assign a transect number to
    %each time point
    adcp(jj).time=datenum(adcp(jj).year+2000,adcp(jj).month,adcp(jj).day,...
        adcp(jj).hour,adcp(jj).minute,zeros(1,length(adcp(jj).minute)));
    datestr(adcp(jj).time(1))
end
