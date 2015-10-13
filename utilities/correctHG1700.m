% this algorithm was implemented by Ifig, 8/20/2013, designed by Zaltan and
% Ifig. Edited by Jianzhu Huai, 8/25/2013. 

% this algorithm estimates the linear mapping between the time interval of
% the raw IMU data and the time interval of the GPS timetagged IMU data.
% Then use this linear model to map the time interval in the raw data to
% get the GPS time increment for each epoch. The time origin was set 
% as the first of the GPS timetag in the GPS timetagged data.

% the basic algorithm: 1, read in the timestamps from the raw HG1700
% records, IMU.txt converted from IMU.bin, and from the GPS timetagged
% HG1700 records in IMUNAV.txt. 2, use linear fit to fit a model between
% the time intervals of the raw data and the gps tagged data, 3, use this
% model to map the raw time interval to the corrected interval. 4, for each
% epoch, its time tag is the corrected accumulated interval plus the time
% origin which is chosen as the first GPS tag. output in IMUCORR.txt.

% the corrected timestamps for HG1700 data improved the performance the
% Kalman filter on HG1700 marginally, but this algorithm takes a lot of
% time
clear all; clc; 
fprintf('Reading data... ');
datadir='C:\Users\huai.3\Desktop\source_testdata_specs\20130719\';
data_1 = dlmread([datadir, 'HG1700IMU.txt'], ',');% this IMU.txt can be obtained by using the readcsv.exe option 4
data_2 = read_hg1700_bin([datadir,'HG1700IMUNAV.bin']); % IMUNAV.bin raw binary file
% data_3= readHG1700navtxt('\\File.ceegs.ohio-state.edu\SPIN\UAV\July 19, 2013\IMU\HG1700\20130719_183532\IMUNAV.txt', 1000);
% data_3 and data_2 are expected to be the same, but data_3 is less
% accurate than data_2.
% data_2=data_2(:,1:16);
% data_1=data_1(1:floor(size(data_1,1)/100),:);
% data_2=data_2(1:floor(size(data_2,1)/100),:);
fprintf('Done\n');

ts_1 = (data_1(:, 16) - data_1(1, 16));
i = 2:size(ts_1, 1);
i2 = i-1;
ts_1_diff = ts_1(i) - ts_1(i2);

ts_2 = data_2(:, 1) - data_2(1, 1);
i = 2:size(ts_2, 1);
i2 = i-1;
ts_2_diff = ts_2(i) - ts_2(i2);

fprintf('Matching timestamps... ');
% Match the timestamps
figure(2); clf; hold on;
pairs = zeros(length(ts_2), 2);
for i = 1 : length(ts_2),
     [minval, ind] = min(abs(ts_1(:)-ts_2(i)));
     pairs(i, 1) = ts_2(i);
     pairs(i, 2) = ts_1(ind);
 end;
 %save pairs;
 
 fprintf('Done\n');
%pnum = 1 : 10: length(pre_ts_1);
pnum = 900 : 1000;

% Figures
figure(1); clf; hold on;
plot(ts_1(2:length(ts_1)), ts_1_diff, 'b-');
xlabel('t [s]');
ylabel('\Deltat [s]');
title('IMU - Original');

figure(2); clf; hold on;
plot(ts_2(2:length(ts_2)), ts_2_diff, 'r-');
xlabel('t [s]');
ylabel('\Deltat [s]');
title('IMUNAV - GPS Time corrected');

% Histograms
figure(3); clf; hold on;
subplot(2,1,1); hold on;
hist(ts_1_diff, 1000);
xlabel('\Deltat [s]')
ylabel('No. [-]')
title('IMU  - Original')
xlim([0 0.05]);

subplot(2,1,2); hold on;
hist(ts_2_diff, 1000);
xlabel('\Deltat [s]')
ylabel('No. [-]')
title('IMUNAV - GPS Time corrected')
xlim([0 0.05]);

%
figure(4); clf; hold on;
subplot(2,2,1); hold on;
plot(pairs(pnum(1),1), 1, 'r*');
plot(pairs(pnum(1),2), 1, 'b*');
plot(pairs(pnum,1), 1, 'r*');
plot(pairs(pnum,2), 1, 'b*');
xlabel('t [s]');
legend('IMUNAV - GPS Time corrected', 'IMU - Original');
title('IMU and IMUNAV')

subplot(2,2,2); hold on;
plot(pairs(:,2), pairs(:,1), 'b*');
axis equal;

% Calculate some statistics
diff = pairs(:,2) - pairs(:,1);
mean_diff = mean(diff);
std_diff = std(diff);
disp('Descriptive statistical values')
disp(['Mean of the differences: ' num2str(mean_diff)]);
disp(['STD of the differences: ' num2str(std_diff)]);
disp(['Median of the differences: ' num2str(median(diff))]);
disp(['Min of the differences: ' num2str(min(diff))]);
disp(['Max of the differences: ' num2str(max(diff))]);

fprintf('removing outliers...');
% Remove outliers
ind_wo_out = find(abs(diff(:) - mean_diff) < 3 * std_diff);
diff_filt = diff(ind_wo_out);
disp('Without outliers')
disp(['Mean of the differences: ' num2str(mean(diff_filt))]);
disp(['STD of the differences: ' num2str(std(diff_filt))]);
disp(['Median of the differences: ' num2str(median(diff_filt))]);
disp(['Min of the differences: ' num2str(min(diff_filt))]);
disp(['Max of the differences: ' num2str(max(diff_filt))]);
fprintf('Done\n');
% Fitting
pre_ts_1 = pairs(ind_wo_out, 2);
pre_ts_2 = pairs(ind_wo_out, 1);
p = polyfit(pre_ts_1,pre_ts_2,1);
figure(4); hold on;
subplot(2,2,2); hold on;
yfit = @(x)  p(1) * x + p(2);
correct_ts = yfit(pre_ts_1);
plot(pre_ts_1,correct_ts, 'r-');
legend('Corresponding ts', 'Linear fitting');
xlabel('IMU - Original');
ylabel('IMUNAV - GPS Time corrected');
title('Fitting')

figure(4); 
subplot(2,2,3); hold on;
plot(pre_ts_1(pnum(1)), 1, 'r*');
plot(pre_ts_2(pnum(1)), 1, 'b*');
plot(correct_ts(pnum(1)), 1, 'g*');
plot(pre_ts_1(pnum), 1, 'r*');
plot(pre_ts_2(pnum), 1, 'b*');
plot(correct_ts(pnum), 1, 'g*');
xlabel('t [s]');
legend('IMUNAV - GPS Time corrected', 'IMU - Original', 'IMUNAV corrected');
title('Corrected IMUNAV timestamps')

res = sqrt(mean((correct_ts - pre_ts_2).^2));
disp(['Residual [s]: ' num2str(res)]);

subplot(2,2,2); hold on;
plot(pairs(:,2), pairs(:,1), 'b-');
axis equal;

% Save IMU correct 
fprintf('Writing data...')
data_2=data_2(ind_wo_out,:);
data_2(:,1) = correct_ts+ data_2(1, 1);
dlmwrite('IMUCORR.txt', data_2, 'precision', '%16.9f');
fprintf('Done\n')
