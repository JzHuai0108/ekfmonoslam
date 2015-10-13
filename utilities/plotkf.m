function f = plotkf(kf, plotDir);
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
%-------------------------------------------------------------------------------
if (nargin<1) return; end;
[nn,mm] = size(kf);
if (nn<1) return; end;
isPrintToFile = false;
if nargin<2 plotDir = []; end;
if ~isempty(plotDir)
    isPrintToFile = true;
    n = length(plotDir);
    if (plotDir(n)~='\')
        plotDir = [plotDir '\'];
    end
    if ~isdir(plotDir)
        mkdir(plotDir);
    end
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
%close all
f(2) = figure;
plot(kf(:,1)-kf(1,1),kf(:,4),'marker','.')
grid
xlabel('Time [s]')
ylabel('m')
title('Height Profile');
if isPrintToFile
    print(f(2), '-dtiff', [plotDir 'ht']);
end
%close all
f(3) = figure;
plot(kf(:,1)-kf(1,1),kf(:,5:7),'marker','.')
grid
xlabel('Time [s]')
ylabel('m/s')
legend('Vx','Vy','Vz')
title('Velocity');
if isPrintToFile
    print(f(3), '-dtiff', [plotDir 'vel']);
end
%close all
% % note: roll is around pi, make to around 0
% kf(:,8) = kf(:,8);
% loc = find(kf(:,8)>pi);
% kf(loc,8) = kf(loc,8);
f(4) = figure;
plot(kf(:,1)-kf(1,1),kf(:,8:10),'marker','.')
grid
xlabel('Time [s]')
ylabel('degree')
legend('Roll','Pitch','Yaw/Heading')
title('Orientation');
if isPrintToFile
    print(f(4), '-dtiff', [plotDir 'att']);
end
%close all
if mm==19
f(5) = figure;
plot(kf(:,1)-kf(1,1),kf(:,11:13),'marker','.')
grid
xlabel('Time [s]')
ylabel('m')
legend('X','Y','Z')
title('Position RMS in ECEF XYZ');
if isPrintToFile
    print(f(5), '-dtiff', [plotDir 'rms_pos']);
end
%close all
f(6) = figure;
plot(kf(:,1)-kf(1,1),kf(:,14:16),'marker','.')
grid
xlabel('Time [s]')
ylabel('m/s')
legend('Vn','Ve','Vd')
title('Velocity RMS in NED');
if isPrintToFile
    print(f(6), '-dtiff', [plotDir 'rms_vel']);
end
%close all
f(7) = figure;
plot(kf(:,1)-kf(1,1),kf(:,17:19)*180/pi,'marker','.')
grid
xlabel('Time [s]')
ylabel('radian')
legend('Roll','Pitch','Yaw/Heading')
title('Orientation RMS in RPY');
if isPrintToFile
    print(f(7), '-dtiff', [plotDir 'rms_att']);
end
end
%close all
