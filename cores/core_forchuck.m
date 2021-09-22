% pb210 info for chuck

%%
clear all,close all,clc

cd C:\GLOVER\output\myanmar\pb210

F=dir('*0.mat');
cname = vertcat(F.name);
cname = cname(:,1:5);

out = [];
for jj=1:length(F)
    cd C:\GLOVER\output\myanmar\pb210

    load(F(jj).name)

tt(:,1) = vertcat(pb210.midinterval);
tt(:,2) = vertcat(pb210.intervalsize);
tt(:,3) = vertcat(pb210.actatcoll_saltcor);
tt(:,4) = vertcat(pb210.salt_cor_mass);

% load mud frac:
cd(['C:\GLOVER\data\myanmar\Pb210\',cname(jj,1:5)])
fid = fopen([cname(jj,1:5),'_NotebookInfo.csv']);
D = textscan(fid,'%*s %*s %f %f %f %f %f %s','Delimiter',',','HeaderLines',4);
fclose(fid);
mud_frac=D{5};

tt(:,5) = mud_frac;

out = [out;tt];
clear tt
end

%%
clear all,close all,clc
cd C:\GLOVER\output\myanmar\grainsize
load('BR_GrainSize_Cores2')


crs = fieldnames(cores);

out1 = [];
out2 = [];

for jj=1:length(crs)
    for kk=1:length(cores.(crs{jj}))
        samp{kk,1} = vertcat(cores.(crs{jj})(kk).name(1:end-4));
    end
    
    tt(:,1) = vertcat(cores.(crs{jj}).median);
    tt(:,2:4) = vertcat(cores.(crs{jj}).frac);
    
    core_num = crs{jj};
    tt(:,5) = str2num(core_num(3:5));
    
    
    
    out1 = [out1;tt];
    out2 = [out2;samp];
    clear tt samp
end




