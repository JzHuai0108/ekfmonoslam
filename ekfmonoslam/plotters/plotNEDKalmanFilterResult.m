function f = plotNEDKalmanFilterResult(kf, plotDir, rotationDim)
%-------------------------------------------------------------------------------
% f = plotkf(kf, isPrintToFile);
% generate plots for Kalman filter estimates
% input:
% 1) kf is the data loaded into memory
% => use kf = readdata(solFileName, 1+9+9) to load binary file format
% => use kf = load(solFileName) to load text file format
% 2) plotDir => to control where to print a tiff file for each plot or not,
% use empty [] to turn off print
% output:
% f => figure handles
% by Dr. Yudan Yi
% modified date on: 04/08/2013
% modified by Jianzhu Huai on 12/17/2021
% 
%-------------------------------------------------------------------------------
if nargin < 3
    rotationDim = 3;  % 3 euler angles rpy, 4 quaternion xyzw.
end
isPrintToFile = false;
if nargin<2 
    plotDir = []; 
end
if isfolder(plotDir)
    isPrintToFile = true;
end
%-------------------------------------------------------------------------------
max_x = max(abs(kf(:,3)));
max_y = max(abs(kf(:,2)));
isNED = false;
if (max_y>pi||max_x>(2.0*pi))
    isNED = true;
end

f(1) = figure;
plot(kf(:,3),kf(:,2),'marker','.')
grid
axis equal
if (isNED)
    xlabel('East [m]')
    ylabel('North[m]')
else
    xlabel('Lontitude [radian]')
    ylabel('Latitude [radian]')
end
title('Ground Track');
if isPrintToFile
    print(f(1), '-dtiff', [plotDir 'ground_track']);
end

f(2) = figure;
plot(kf(:,1)-kf(1,1),kf(:,4),'marker','.')
grid
xlabel('Time [s]')
ylabel('m')
title('Height Profile');
if isPrintToFile
    print(f(2), '-dtiff', [plotDir 'ht']);
end

f(3) = figure;
plot(kf(:,1)-kf(1,1),kf(:,5:7),'marker','.');
grid
xlabel('Time [s]')
ylabel('m/s')
legend('Vx','Vy','Vz')
title('Velocity');
if isPrintToFile
    print(f(3), '-dtiff', [plotDir 'vel']);
end

% % note: roll is around pi, make to around 0
% kf(:,8) = kf(:,8);
% loc = find(kf(:,8)>pi);
% kf(loc,8) = kf(loc,8);
nextFig=4;
nextDim=8;
f(nextFig) = figure;
if rotationDim == 3
    plot(kf(:,1)-kf(1,1),kf(:, nextDim+(0:2)),'marker','.');
else
    eulerangles = zeros(size(kf, 1), 3);
    for i = 1 : size(kf, 1)
        N_q_S = [kf(i, nextDim + 3), kf(i, nextDim + (0:2))]';
        eulerradians = rotqr2eu('xyz', N_q_S);
        eulerangles(i, :) = eulerradians' *180/pi;
    end
    plot(kf(:,1)-kf(1,1), eulerangles, 'marker','.');
end
grid;
xlabel('Time [s]');
ylabel('degree');
legend('Roll','Pitch','Yaw/Heading');
title('Orientation');
if isPrintToFile
    print(f(nextFig), '-dtiff', [plotDir 'att']);
end

nextFig=nextFig+1;
nextDim=nextDim+rotationDim;
f(nextFig) = figure;
plot(kf(:,1)-kf(1,1),kf(:,nextDim+(0:2)),'marker','.');
grid;
xlabel('Time [s]');
ylabel('m');
legend('N','E','D');
title('Position RMS in NED');
if isPrintToFile
    print(f(nextFig), '-dtiff', [plotDir 'rms_pos']);
end

nextFig=nextFig+1;
nextDim=nextDim+3;
f(nextFig) = figure;
plot(kf(:,1)-kf(1,1),kf(:,nextDim+(0:2)),'marker','.');
grid;
xlabel('Time [s]');
ylabel('m/s');
legend('Vn','Ve','Vd');
title('Velocity RMS in NED');
if isPrintToFile
    print(f(nextFig), '-dtiff', [plotDir 'rms_vel']);
end

nextFig=nextFig+1;
nextDim=nextDim+3;
f(nextFig) = figure;
plot(kf(:,1)-kf(1,1),kf(:,nextDim+(0:2))*180/pi,'marker','.');
grid;
xlabel('Time [s]');
ylabel('degree');
legend('Roll','Pitch','Yaw/Heading');
title('Orientation RMS in RPY');
if isPrintToFile
    print(f(nextFig), '-dtiff', [plotDir 'rms_att']);
end

