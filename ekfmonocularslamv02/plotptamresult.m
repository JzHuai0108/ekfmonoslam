function plotptamresult()
clc
addpath('..\voicebox\'); % for rotro2qr and rotqr2eu, they are more robuts
fn='H:\relaylatest\EKF_monoSLAM_1pRANSAC\sequences\ic\DeletemeTrans.txt';
alldata=load(fn);
figure(2)
for jerk=40:size(alldata);
plot3( alldata(1:jerk,5),alldata(1:jerk,6),alldata(1:jerk,7), '-k');
hold on
row2c=roteu2ro('xyz',alldata(jerk, 2:4));
qc2w=rotro2qr(row2c');
draw_camera( [alldata(jerk,5:7)'; qc2w], 'k');
hold off
title('Tcinw');
xlabel('x');
ylabel('y');
zlabel('z');
legend('x', 'y', 'z');
view(0,0);
pause(0.05);
end