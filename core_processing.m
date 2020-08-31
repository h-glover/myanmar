%%
clear all,close all,clc

cname = 'BR432';
cd(['C:\GLOVER\data\myanmar\Pb210\',cname])

% input the Supported Activity (start with 0 and examine data), should be
% less than 1 and constant for all profiles in an area
supp_act=0.9;

% First use PeakCounter to output the raw counts/data from the SPE files
cts=PeakCounter('.',0);

% load the sampling data from the excel file
fid = fopen([cname,'_NotebookInfo.csv']);
D = textscan(fid,'%*s %*s %f %f %f %f %f %s','Delimiter',',','HeaderLines',4);fclose(fid);
tin=D{1};
wet=D{2};
dry=D{3};
weight=D{4};
mud_frac=D{5};
plate_date=datenum(D{6});
% get the sample collection date from the file as well
fid = fopen([cname,'_NotebookInfo.csv']);
coll_date=textscan(fid, '%*s %s[^\n]', 1, 'HeaderLines',1,'Delimiter',',');fclose(fid);
coll_date=datenum(coll_date{1});

% process all of the depths
for jj=1:length(cts)
    pb210(jj) = ...
        pb210_proc(cts(jj),plate_date(jj),coll_date,weight(jj),...
        supp_act,mud_frac(jj),tin(jj),wet(jj),dry(jj));
end


% save the output for later use and plotting:
output_name = ['C:\GLOVER\output\myanmar\pb210\',cname,'_pb210ex'];
save(output_name,'pb210')