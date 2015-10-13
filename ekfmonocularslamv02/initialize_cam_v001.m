%-----------------------------------------------------------------------
% 1-point RANSAC EKF SLAM from a monocular sequence
%-----------------------------------------------------------------------

% Copyright (C) 2010 Javier Civera and J. M. M. Montiel
% Universidad de Zaragoza, Zaragoza, Spain.

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation. Read http://www.gnu.org/copyleft/gpl.html for details

% If you use this code for academic work, please reference:
%   Javier Civera, Oscar G. Grasa, Andrew J. Davison, J. M. M. Montiel,
%   1-Point RANSAC for EKF Filtering: Application to Real-Time Structure from Motion and Visual Odometry,
%   to appear in Journal of Field Robotics, October 2010.

%-----------------------------------------------------------------------
% Authors:  Javier Civera -- jcivera@unizar.es
%           J. M. M. Montiel -- josemari@unizar.es

% Robotics, Perception and Real Time Group
% Aragón Institute of Engineering Research (I3A)
% Universidad de Zaragoza, 50018, Zaragoza, Spain
% Date   :  May 2010
%-----------------------------------------------------------------------

function [cam, resizeScale] = initialize_cam(camtype, maxEdge)
if(nargin<2)
    maxEdge=1000;
end
resizeScale=0;
switch camtype
    case 0
        % Civera's classic model as in Civera, 2008
        d =     0.0112;
        cam.nRows = 240;
        cam.nCols = 320;
        Cx =    1.7945 / d;
        Cy =    1.4433 / d;
        k1=     6.333e-2;
        k2=     1.390e-2;
        f =     2.1735;     
        % The equivalent calibration parameters in Bouguet's format
        % provided also by Civera's mat file
        cam.fc=[ 195.547279088826;195.547279088826 ];
        cam.cc=[173.626763984641;131.116786110196]; 
        cam.kc=[-0.342121122410104;0.146843863769195;0.000850015447411755;0;-0.0319609120043733];
        cam.alpha=0;        
        cam.K =     sparse( [ cam.fc(1)   0     cam.cc(1);
            0  cam.fc(2)    cam.cc(2);
            0    0     1] );
        
        cam.model = 'two_distortion_parameters';
    case 1% casio 2 on experiment aug 08 2013 on GPSvan
        cam.nRows = 720;
        cam.nCols = 1280;
        edgeLeg=max(cam.nRows, cam.nCols);
        while(edgeLeg>maxEdge)
            edgeLeg=round(edgeLeg/2);
            resizeScale=resizeScale+1;
        end
        scale=2^resizeScale;
        cam.nRows = 720/scale;
        cam.nCols = 1280/scale;
        
        cam.fc=[975.926;976.2174]/scale;     
        cam.cc=[646.4251;368.0766]/scale; 
        cam.kc=[-0.17483	0.20382	0.00189	0.00329 0]';
        cam.alpha=0;     
        
        cam.K =     sparse( [ cam.fc(1)   cam.alpha*cam.fc(1)   cam.cc(1);
            0  cam.fc(2)    cam.cc(2);
            0    0     1] );        
        cam.model = 'five_distortion_parameters';
    case 2% nikon N800
        d=0.00488; %pixel size in x and y directionis    
        cam.nRows = 4912;
        cam.nCols = 7360;
        edgeLeg=max(cam.nRows, cam.nCols);
        while(edgeLeg>maxEdge)
            edgeLeg=round(edgeLeg/2);
            resizeScale=resizeScale+1;
        end
        scale=2^resizeScale;        
        cam.nRows = 4912/scale;
        cam.nCols = 7360/scale;
        Cx =    3733.85248/scale;
        Cy =   2387.86017/scale;
        k1=     -0.19954;
        k2=     0.23275;
        p1=-0.00021;
        p2=0.00060;% notation in accord with Jean-Yves Bouguet
       
        fx=11021.99028/scale;
        fy=11026.76793/scale;
        cam.fc=[fx;fy];     
        cam.cc=[Cx;Cy]; 
        cam.alpha=0;    
        cam.kc=[k1,k2,p1,p2,0]';
        cam.K =     sparse( [ fx   cam.alpha*fx   Cx;
            0  fy    Cy;
            0    0     1] );        
        cam.model = 'five_distortion_parameters';   
    case 3% casio 2 on experiment aug 08 2013 on GPSvan disabling the distortion
        cam.nRows = 720;
        cam.nCols = 1280;
        edgeLeg=max(cam.nRows, cam.nCols);
        while(edgeLeg>maxEdge)
            edgeLeg=round(edgeLeg/2);
            resizeScale=resizeScale+1;
        end
        scale=2^resizeScale;
        cam.nRows = 720/scale;
        cam.nCols = 1280/scale;
        
        cam.fc=[975.926;976.2174]/scale;     
        cam.cc=[646.4251;368.0766]/scale; 
        cam.kc=zeros(5,1); %[-0.17483	0.20382	0.00189	0.00329 0]';
        cam.alpha=0;           
        cam.K =     sparse( [ cam.fc(1)   cam.alpha*cam.fc(1)   cam.cc(1);
            0  cam.fc(2)    cam.cc(2);
            0    0     1] );        
        cam.model = 'five_distortion_parameters';
    otherwise
end
% cam.k1 =    k1;
% cam.k2 =    k2;
% cam.Cx =    Cx;
% cam.Cy =    Cy;
% cam.f =     f;
% cam.dx =    d;
% cam.dy =    d;