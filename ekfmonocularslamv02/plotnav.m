% plot the filtered estimates for navigation states, imu errors and imu
% error covariance stds.

kf = readdata(filresfile, 1+18);
if (~isOutNED)
%    save_google_kml(kf(:,2:4), kmlfilename);
end
if(useGPS)
inillh_ant2=options.inillh_ant;

posdata=loadAllGPSData(gpsfile, [kf(1,1), kf(end,1)], gpspostype);
maxl=size(posdata,1);
if(isOutNED)
    for i=1:maxl
        rovXYZ = posdata(i,2:4);
        posdata(i,2:4)=posdiff_v001(rovXYZ',inillh_ant2);
    end
end
end
nextFig=1;
f(nextFig) = figure;
if (isOutNED)
    %         plot(kf(:,3),kf(:,2),'g.')
    %         hold on
    %         plot(posdata(:,3),posdata(:,2),'r+')
    %         grid
    %         xlabel('East [m]')
    %         ylabel('North[m]')
    plot3(kf(:,1)-kf(1,1),kf(:,3),kf(:,2),'g.')
    hold on
    if(useGPS)
    plot3(posdata(:,1)-kf(1,1),posdata(:,3),posdata(:,2),'r+')
    end
    grid
    axis equal
    xlabel('Time [s]')
    ylabel('East [m]')
    zlabel('North[m]')
else
    plot(kf(:,3),kf(:,2),'g.')
    hold on
    if(useGPS)
    plot(posdata(:,3),posdata(:,2),'r+')
    end
    grid
    axis equal
    xlabel('Longitude [radian]')
    ylabel('Lattitude [radian]')
end
title('green KF trajectory and the red GPS reference for GPS antenna');
saveas(f(nextFig),[resdir,'red truth and track'],'fig');

nextFig=nextFig+1;
f(nextFig) = figure;
plot(kf(:,1)-kf(1,1),kf(:,4),'g.')
hold on
if(useGPS)
plot(posdata(:,1)-kf(1,1),posdata(:,4),'r+');
end
grid
xlabel('Time [s]')
ylabel('Height/m')
title('Height of antenna by KF(green) and reference(red)');
saveas(f(nextFig),[resdir 'red truth and height'],'fig');

isPrintToFile=0;
if(options.useCam)
    campose=load(options.camPoseFile);
    if(~isempty(campose))
        nextFig=nextFig+1;
        f(nextFig) = figure;
        plot(campose(:,2)-campose(1,2),campose(:,6:8),'marker','.')
        grid
        xlabel('Time [s]')
        ylabel('m')
        legend('x','y','z')
        title('Ts2c');
        if isPrintToFile
            print(f(nextFig), '-dtiff', [plotDir 'Ts2c']);
        end
        
        nextFig=nextFig+1;
        f(nextFig) = figure;
        plot(campose(:,2)-campose(1,2),campose(:, 3:5),'marker','+');
        grid
        xlabel('Time [s]')
        ylabel('degree')
        legend('Roll','Pitch','Yaw/Heading')
        title('Cs2c');
        if isPrintToFile
            print(f(nextFig), '-dtiff', [plotDir 'qs2c']);
        end
        nextFig=nextFig+1;
        f(nextFig) = figure;
        plot(campose(:,2)-campose(1,2),campose(:,12:14),'marker','.')
        grid
        xlabel('Time [s]')
        ylabel('m')
        legend('x','y','z')
        title('Xs02e');
        if isPrintToFile
            print(f(nextFig), '-dtiff', [plotDir 'Xs02e']);
        end
        
        nextFig=nextFig+1;
        f(nextFig) = figure;
        plot(campose(:,2)-campose(1,2),campose(:, 9:11),'marker','+');
        grid
        xlabel('Time [s]')
        ylabel('degree')
        legend('Roll','Pitch','Yaw/Heading')
        title('qs02e');
        if isPrintToFile
            print(f(nextFig), '-dtiff', [plotDir 'qs02e']);
        end
    end
end
%close all
filresfile=[resdir, 'filresult.bin'];
kf = readdata(filresfile, 1+18);
plotkf_v001(kf, resdir);
imuresfile=[resdir, 'imuresult.bin'];
err = readdata(imuresfile, 1+24);
ploterr_v001(err,resdir);