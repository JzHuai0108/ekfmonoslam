clear;
load input.mat; %read the input
load input_adi.mat;
fs=10;

%compute the Allan variance curve
[t_axis, allan, Nrw]=cp_allan_var_v000(input,fs,1,60);
[t_axis_adi, allan_adi, Nrw_adi]=cp_allan_var_v000(input_adi,fs,2,60);

%approximation and rw parameter
apx_lvl=16;
[inp_appx, boun_appx]=lvl_filter_v000(input, apx_lvl,0,0);
rw_coef=cp_wn_pow_v000(input,apx_lvl+1,0);
rw_pow=diag(rw_coef'*rw_coef).^0.5;

apx_lvl_adi=16;
[inp_appx_adi, boun_appx_adi]=lvl_filter_v000(input_adi, apx_lvl_adi,0,0);
rw_coef_adi=cp_wn_pow_v000(input_adi,apx_lvl_adi+1,0);
rw_pow_adi=diag(rw_coef_adi'*rw_coef_adi).^0.5;

%White noise level
whi_coef=cp_wn_pow_v000(input,3,1);
whi_pow=diag(whi_coef'*whi_coef).^0.5;
whi_coef_adi=cp_wn_pow_v000(input_adi,3,1);
whi_pow_adi=diag(whi_coef_adi'*whi_coef_adi).^0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Individual IMU Results %%%%%%%%%%
%%%%% AMS
flt_lvl=9;
[inp_filt, boun_filt]=lvl_filter_v000(input,flt_lvl,0,0);
inputx=inp_filt-inp_appx;
cor_sz=2^16;
inp_corr_x=xcorr(inputx(1,:),cor_sz,'biased');
inp_corr_y=xcorr(inputx(2,:),cor_sz,'biased');
inp_corr_z=xcorr(inputx(3,:),cor_sz,'biased');

plot(0:cor_sz,inp_corr_x(cor_sz+1:2*cor_sz+1));grid;
hold on;
plot(0:cor_sz,inp_corr_y(cor_sz+1:2*cor_sz+1),'r');
plot(0:cor_sz,inp_corr_z(cor_sz+1:2*cor_sz+1),'k');

inp_corr=[inp_corr_x(:,cor_sz+1:2*cor_sz+1);inp_corr_y(:,cor_sz+1:2*cor_sz+1);inp_corr_z(:,cor_sz+1:2*cor_sz+1)];
[mod_h exp_h]=expo_fit_v000([(2^flt_lvl),8000;(2^flt_lvl),15000;(2^flt_lvl),20000], inp_corr,1,1);

flt_lvl=4;
[inp_filt, boun_filt]=lvl_filter_v000(input,flt_lvl,0,0);
inputx=inp_filt-inp_appx;
cor_sz=2^11;
inp_corr_x=xcorr(inputx(1,:),cor_sz,'biased');
inp_corr_y=xcorr(inputx(2,:),cor_sz,'biased');
inp_corr_z=xcorr(inputx(3,:),cor_sz,'biased');

prev_cor_x=exp_h(1,2)*exp_h(1,1).^[0:cor_sz]+exp_h(1,3);
prev_cor_y=exp_h(2,2)*exp_h(2,1).^[0:cor_sz]+exp_h(2,3);
prev_cor_z=exp_h(3,2)*exp_h(3,1).^[0:cor_sz]+exp_h(3,3);
figure;
plot(0:cor_sz,inp_corr_x(cor_sz+1:2*cor_sz+1)-prev_cor_x);grid;
hold on;
plot(0:cor_sz,inp_corr_y(cor_sz+1:2*cor_sz+1)-prev_cor_y,'r');
plot(0:cor_sz,inp_corr_z(cor_sz+1:2*cor_sz+1)-prev_cor_z,'k');

inp_corr=[inp_corr_x(:,cor_sz+1:2*cor_sz+1)-prev_cor_x;inp_corr_y(:,cor_sz+1:2*cor_sz+1)-prev_cor_y;inp_corr_z(:,cor_sz+1:2*cor_sz+1)-prev_cor_z];
[mod_l exp_l]=expo_fit_v000([16,500;16, 1200;16, 1200], inp_corr,1,1);

%compare the allan results
allan_th=arma2allan_v000([1 -mod_h(1,1);1 -mod_l(1,1);1 0;1 -1],[mod_h(1,2);mod_l(1,2);whi_pow(1);rw_pow(1)],fs,t_axis,allan(1,:));
allan_th=arma2allan_v000([1 -mod_h(2,1);1 -mod_l(2,1);1 0;1 -1],[mod_h(2,2);mod_l(2,2);whi_pow(2);rw_pow(2)],fs,t_axis,allan(2,:));
allan_th=arma2allan_v000([1 -mod_h(3,1);1 -mod_l(3,1);1 0;1 -1],[mod_h(3,2);mod_l(3,2);whi_pow(3);rw_pow(3)],fs,t_axis,allan(3,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ADI
flt_lvl=9;
[inp_filt, boun_filt]=lvl_filter_v000(input_adi,flt_lvl,0,0);
inputx=inp_filt-inp_appx_adi;
cor_sz=2^16;
inp_corr_x=xcorr(inputx(1,:),cor_sz,'biased');
inp_corr_y=xcorr(inputx(2,:),cor_sz,'biased');
inp_corr_z=xcorr(inputx(3,:),cor_sz,'biased');

figure;
plot(0:cor_sz,inp_corr_x(cor_sz+1:2*cor_sz+1));grid;
hold on;
plot(0:cor_sz,inp_corr_y(cor_sz+1:2*cor_sz+1),'r');
plot(0:cor_sz,inp_corr_z(cor_sz+1:2*cor_sz+1),'k');

inp_corr=[inp_corr_x(:,cor_sz+1:2*cor_sz+1);inp_corr_y(:,cor_sz+1:2*cor_sz+1);inp_corr_z(:,cor_sz+1:2*cor_sz+1)];
[mod_h_adi exp_h_adi]=expo_fit_v000([(2^flt_lvl),18000;(2^flt_lvl),18000;(2^flt_lvl),18000], inp_corr,1,1);

%compare the allan results
allan_th=arma2allan_v000([1 -mod_h_adi(1,1);1 0;1 -1],[mod_h_adi(1,2);whi_pow_adi(1);rw_pow_adi(1)],fs,t_axis,allan_adi(1,:));
allan_th=arma2allan_v000([1 -mod_h_adi(2,1);1 0;1 -1],[mod_h_adi(2,2);whi_pow_adi(2);rw_pow_adi(2)],fs,t_axis,allan_adi(2,:));
allan_th=arma2allan_v000([1 -mod_h_adi(3,1);1 0;1 -1],[mod_h_adi(3,2);whi_pow_adi(3);rw_pow_adi(3)],fs,t_axis,allan_adi(3,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Combined (Multi) IMU Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% combine outputs
R=diag([diag(whi_coef'*whi_coef);diag(whi_coef_adi'*whi_coef_adi)]);
H=[sin(45*pi/180) cos(45*pi/180);
   sin(90*pi/180) cos(90*pi/180);
   sin(70*pi/180) cos(70*pi/180)
   sin(10*pi/180) cos(10*pi/180)
   sin(0*pi/180) cos(0*pi/180)
   sin(-10*pi/180) cos(-10*pi/180)];
%standard method
T2=inv(H'*inv(R)*H)*H'*inv(R);
output=T2*[input;input_adi];

%compute the Allan variance curve
[t_axis, allan_out, Nrw]=cp_allan_var_v000(output,fs,3,60);

%approximation and rw parameter
apx_lvl=16;
[out_appx, boun_appx]=lvl_filter_v000(output, apx_lvl,0,0);
rw_coef_out=cp_wn_pow_v000(output,apx_lvl+1,0);
rw_pow_out=diag(rw_coef_out'*rw_coef_out).^0.5;

%White noise level
whi_coef_out=cp_wn_pow_v000(output,3,1);
whi_pow_out=diag(whi_coef_out'*whi_coef_out).^0.5;

%%%%%%% Start modeling
flt_lvl=9;
llim=2^flt_lvl; %autocorr is valid for k>=llim;
[out_filt, boun_filt]=lvl_filter_v000(output,flt_lvl,0,0);
outputx=out_filt-out_appx;
cor_sz=2^15;
out_corr_x=xcorr(outputx(1,:),cor_sz,'biased');
out_corr_y=xcorr(outputx(2,:),cor_sz,'biased');
out_corr_xy=xcorr(outputx(1,:),outputx(2,:),cor_sz,'biased');

figure;
plot(0:cor_sz,out_corr_x(cor_sz+1:2*cor_sz+1));grid;
hold on;
plot(0:cor_sz,out_corr_y(cor_sz+1:2*cor_sz+1),'r');
plot(0:cor_sz,out_corr_xy(cor_sz+1:2*cor_sz+1),'k');
plot(-cor_sz:0,out_corr_xy(1:cor_sz+1),'k');

out_corr=[out_corr_x(:,cor_sz+1:2*cor_sz+1);out_corr_y(:,cor_sz+1:2*cor_sz+1);out_corr_xy(:,cor_sz+1:2*cor_sz+1);fliplr(out_corr_xy(:,1:cor_sz+1));];
[mod_h_out exp_h_out]=expo_fit_v000([ones(3,1)*[llim,20000];[llim 10000]], out_corr,1,1);

%%%lower level
flt_lvl=3;
[out_filt, boun_filt]=lvl_filter_v000(output,flt_lvl,0,0);
outputx=out_filt-out_appx;
cor_sz=2^10;
out_corr_x=xcorr(outputx(1,:),cor_sz,'biased');
out_corr_y=xcorr(outputx(2,:),cor_sz,'biased');
out_corr_xy=xcorr(outputx(1,:),outputx(2,:),cor_sz,'biased');

prev_cor_x=exp_h_out(1,2)*exp_h_out(1,1).^[0:cor_sz]+exp_h_out(1,3);
prev_cor_y=exp_h_out(2,2)*exp_h_out(2,1).^[0:cor_sz]+exp_h_out(2,3);
prev_cor_xy_pos=exp_h_out(3,2)*exp_h_out(3,1).^[0:cor_sz]+exp_h_out(3,3);
prev_cor_xy_neg=exp_h_out(4,2)*exp_h_out(4,1).^abs([-cor_sz:0])+exp_h_out(4,3);

plot(0:cor_sz,out_corr_x(cor_sz+1:2*cor_sz+1)-prev_cor_x);grid;
hold on;
plot(0:cor_sz,out_corr_y(cor_sz+1:2*cor_sz+1)-prev_cor_y,'r');
plot(0:cor_sz,out_corr_xy(cor_sz+1:2*cor_sz+1)-prev_cor_xy_pos,'k');
plot(-cor_sz:0,out_corr_xy(1:cor_sz+1)-prev_cor_xy_neg,'k');

out_corr=[out_corr_x(:,cor_sz+1:2*cor_sz+1)-prev_cor_x;out_corr_y(:,cor_sz+1:2*cor_sz+1)-prev_cor_y;out_corr_xy(:,cor_sz+1:2*cor_sz+1)-prev_cor_xy_pos;fliplr(out_corr_xy(:,1:cor_sz+1)-prev_cor_xy_neg)];
[mod_l_out exp_l_out]=expo_fit_v000(ones(4,1)*[8,1000], out_corr,1,1);

allan_th=arma2allan_v000([1 -mod_h_out(1,1);1 -mod_l_out(1,1);1 0;1 -1],[mod_h_out(1,2);mod_l_out(1,2);whi_pow_out(1);rw_pow_out(1)],fs,t_axis,allan_out(1,:));
allan_th=arma2allan_v000([1 -mod_h_out(2,1);1 -mod_l_out(2,1);1 0;1 -1],[mod_h_out(2,2);mod_l_out(2,2);whi_pow_out(2);rw_pow_out(2)],fs,t_axis,allan_out(2,:));


%%%Readjust parameters for multi model
%% 1st order MW appx
[exp_h_apx]=exp_approx_v000(ones(4,1)*[2^9,15000], exp_h_out,1);
[exp_l_apx]=exp_approx_v000(ones(4,1)*[1,1000], exp_l_out,1);

mod_h_apx=exp2markov_v000(exp_h_apx(1:3,:),1);
mod_l_apx=exp2markov_v000(exp_l_apx(1:3,:),1);

allan_th=arma2allan_v000([1 -mod_h_apx(1,1);1 -mod_l_apx(1,1);1 0;1 -1],[mod_h_apx(1,2);mod_l_apx(1,2);whi_pow_out(1);rw_pow_out(1)],fs,t_axis,allan_out(1,:));
allan_th=arma2allan_v000([1 -mod_h_apx(2,1);1 -mod_l_apx(2,1);1 0;1 -1],[mod_h_apx(2,2);mod_l_apx(2,2);whi_pow_out(2);rw_pow_out(2)],fs,t_axis,allan_out(2,:));

%% 2nd method (mod_l will be approximated as second order process)
exp_h2_apx=exp_approx_v000(ones(4,1)*[2^9,15000], exp_h_out,1);
exp_l2a_apx=exp_approx_v000(ones(2,1)*[10,1000], exp_l_out([1 3],:),1);
exp_l2b_apx=[exp_l2a_apx(1,1) (exp_l2a_apx(2,2)^2)/exp_l2a_apx(1,2)];
[mod_l2bres_apx exp_l2bres_apx]=exp_approx2_v000([0,1000], exp_l_out(2,:), exp_l2b_apx,1);

[mod_h2_apx mx_a]=exp2markov_v000(exp_h2_apx(1:3,:),1);
[mod_l2a_apx mx_a]=exp2markov_v000([exp_l2a_apx(1,:);exp_l2b_apx;exp_l2a_apx(2,:)],1);

allan_th=arma2allan_v000([1 -mod_h2_apx(1,1);1 -mod_l2a_apx(1,1);1 0;1 -1],[mod_h2_apx(1,2);mod_l2a_apx(1,2);whi_pow_out(1);rw_pow_out(1)],fs,t_axis,allan_out(1,:));
allan_th=arma2allan_v000([1 -mod_h2_apx(2,1) 0;1 -mod_l2a_apx(2,1) 0; 1 -mod_l2bres_apx(1,1) -mod_l2bres_apx(1,2);1 0 0;1 -1 0],[mod_h2_apx(2,2);mod_l2a_apx(2,2);mod_l2bres_apx(1,3);whi_pow_out(2);rw_pow_out(2)],fs,t_axis,allan_out(2,:));
