function f = ploterr_v001(err, plotDir)
%-------------------------------------------------------------------------------
% f = ploterr(err, isPrintToFile);
% generate plots for Kalman filter error estimates
% input: 
% 1) err is the data loaded into memory
% => use err = readdata(errFileName, 1+15+15) to load binary file format
% => use err = load(errFileName) to load text file format
% 2) plotDir => to control where to print a tiff file for each plot or not,
% use empty [] to turn off print
% output:
% f => figure handles
% by Dr. Yudan Yi
% modified date on: 04/08/2013
%-------------------------------------------------------------------------------
if (nargin<1) return; end;
[nn,mm] = size(err);
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
f(1) = figure;
loc = 2:4;
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
ylabel('m/s/s')
legend('X','Y','Z')
title('Accelerometer bias drift');
if isPrintToFile
    print(f(1), '-dtiff', [plotDir 'bda']);
end
%close all
f(2) = figure;
loc = 5:7;
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
ylabel('radian/s')
legend('X','Y','Z')
title('Gyro bias drift');
if isPrintToFile
    print(f(2), '-dtiff', [plotDir 'bdg']);
end
%close all

f(3) = figure;
loc = 8:10;
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
legend('X','Y','Z')
title('Accelerometer scale factor error');
if isPrintToFile
    print(f(3), '-dtiff', [plotDir 'sfa']);
end
%close all
f(4) = figure;
loc = 11:13;
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
legend('X','Y','Z')
title('Gyro scale factor error');
if isPrintToFile
    print(f(4), '-dtiff', [plotDir 'sfg']);
end
%close all
nextFig=5;
nextDim=14;
%close all
f(nextFig) = figure;
loc = nextDim+(0:2);
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
legend('X','Y','Z')
title('Accelerometer bias drift cov std');
if isPrintToFile
    print(f(nextFig), '-dtiff', [plotDir 'rms_bda']);
end
%close all
nextFig=nextFig+1;
nextDim=nextDim+3;
f(nextFig) = figure;
loc = nextDim+(0:2);
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
legend('X','Y','Z')
title('Gyro bias drift cov std');
if isPrintToFile
    print(f(nextFig), '-dtiff', [plotDir 'rms_bdg']);
end
nextFig=nextFig+1;
nextDim=nextDim+3;
f(nextFig) = figure;
loc = nextDim+(0:2);
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
legend('X','Y','Z')
title('Accelerometer scale factor cov std');
if isPrintToFile
    print(f(nextFig), '-dtiff', [plotDir 'rms_sfa']);
end
%close all
nextFig=nextFig+1;
nextDim=nextDim+3;
f(nextFig) = figure;
loc = nextDim+(0:2);
plot(err(:,1)-err(1,1),err(:,loc),'marker','.')
grid
xlabel('Time [s]')
legend('X','Y','Z')
title('Gyro scale factor error cov std');
if isPrintToFile
    print(f(nextFig), '-dtiff', [plotDir 'rms_sfg']);
end
%close all