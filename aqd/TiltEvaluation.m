% Nan bad data when patch+roll>30
% Fix pressure
% clean mat file by finding when instrument is in water
clear all,close all
load('MMsm_Sept17_aqd.mat')
figure;
subplot(411)
plot(aqd.head),ylabel('heading'),ylim([0 360]),title('Small Sept17')
subplot(412)
plot(aqd.pitch),ylabel('pitch')
subplot(413)
plot(aqd.roll),ylabel('roll')
subplot(414)
plot(abs(aqd.pitch)+abs(aqd.roll)),ylabel('abs pitch+roll')
r=refline(0,10);r.Color='k';r=refline(0,30);r.Color='k';

load('MMLG_Sept17_aqd.mat')
figure;
subplot(411)
plot(aqd.head),ylabel('heading'),ylim([0 360]),title('Large Sept17')
subplot(412)
plot(aqd.pitch),ylabel('pitch')
subplot(413)
plot(aqd.roll),ylabel('roll')
subplot(414)
plot(abs(aqd.pitch)+abs(aqd.roll)),ylabel('abs pitch+roll')
r=refline(0,10);r.Color='k';r=refline(0,30);r.Color='k';

load('MMsm_Mar18_aqd.mat')
figure;
subplot(411)
plot(aqd.head),ylabel('heading'),ylim([0 360]),title('Small Mar18')
subplot(412)
plot(aqd.pitch),ylabel('pitch')
subplot(413)
plot(aqd.roll),ylabel('roll')
subplot(414)
plot(abs(aqd.pitch)+abs(aqd.roll)),ylabel('abs pitch+roll')
r=refline(0,10);r.Color='k';r=refline(0,30);r.Color='k';


load('MMLG_Mar18_aqd.mat')
figure;
subplot(411)
plot(aqd.head),ylabel('heading'),ylim([0 360]),title('Large Mar18')
subplot(412)
plot(aqd.pitch),ylabel('pitch')
subplot(413)
plot(aqd.roll),ylabel('roll')
subplot(414)
plot(abs(aqd.pitch)+abs(aqd.roll)),ylabel('abs pitch+roll')
r=refline(0,10);r.Color='k';r=refline(0,30);r.Color='k';

load('MMAG_Mar18_aqd.mat')
figure;
subplot(411)
plot(aqd.head),ylabel('heading'),ylim([0 360]),title('Agri Mar18')
subplot(412)
plot(aqd.pitch),ylabel('pitch')
subplot(413)
plot(aqd.roll),ylabel('roll')
subplot(414)
plot(abs(aqd.pitch)+abs(aqd.roll)),ylabel('abs pitch+roll')
r=refline(0,10);r.Color='k';r=refline(0,30);r.Color='k';
export_fig AgriChannel_Mar18