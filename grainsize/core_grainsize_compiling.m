clear all,close all,clc

cd C:\GLOVER\data\myanmar\GrainSize_Processing\cores


F=dir('*.csv');
D=[];
for jj=1:length(F)
    temp = ls13320_read(F(jj).name);
    idx(jj) = str2num(F(jj).name(3:5));
    D=[D;temp];temp=[];
end
crs = unique(idx);
cores.br323 = D(idx==crs(1));
cores.br324 = D([idx==crs(2)]);
cores.br326 = D([idx==crs(3)]);
cores.br335 = D([idx==crs(4) | idx==crs(5)]);
cores.br360 = D([idx==crs(6)]);
cores.br430 = D([idx==crs(7)]);
cores.br432 = D([idx==crs(8)]);


cd C:\GLOVER\output\myanmar\grainsize
save('BR_GrainSize_Cores','cores')
