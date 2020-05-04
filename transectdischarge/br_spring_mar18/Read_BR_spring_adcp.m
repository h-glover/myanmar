%% Read in ADCP data (Legacy ASCII text files from WinRiver).  Requires that you have the function adcp_ascii (D. Nowacki)

% Folder containing only text files of each transect
FileList=dir('*.TXT');




for n = 1:length(FileList)
    
    adcp(n) = adcp_ascii(FileList(n).name);

    disp(['processing ' num2str(n)]);
end


% if you're satisfied with the matric produced, make sure to save it in
% your workspace!