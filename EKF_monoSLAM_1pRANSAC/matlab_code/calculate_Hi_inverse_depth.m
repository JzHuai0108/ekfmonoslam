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

function Hi = calculate_Hi_inverse_depth( Xv_km1_k, yi, camera, i, features_info )

zi = features_info(i).h;

number_of_features = size( features_info, 2 );
inverse_depth_features_index = zeros( number_of_features, 1 );
cartesian_features_index = zeros( number_of_features, 1 );

for j=1:number_of_features
    if strncmp(features_info(j).type, 'inversedepth', 1)
        inverse_depth_features_index(j) = 1;
    end
    if strncmp(features_info(j).type, 'cartesian', 1)
        cartesian_features_index(j) = 1;
    end
end

Hi = zeros(2, 13+3*sum(cartesian_features_index)+6*sum(inverse_depth_features_index));

Hi(:,1:13) = dh_dxv( camera, Xv_km1_k, yi, zi );

index_of_insertion = 13 + 3*sum(cartesian_features_index(1:i-1)) + 6*sum(inverse_depth_features_index(1:i-1))+1;
Hi(:,index_of_insertion:index_of_insertion+5) = dh_dy( camera, Xv_km1_k, yi, zi );

return



function Hii = dh_dy( camera, Xv_km1_k, yi, zi )

    Hii = dh_dhrl( camera, Xv_km1_k, yi, zi ) * dhrl_dy( Xv_km1_k, yi );
    
%     % Verification, OK
%         [X,FVAL,EXITFLAG,OUTPUT,J] = fsolve( 'F_Test_dh_dy_newPar', yi, ...
%         optimset( 'Display', 'off', ...
%         'NonlEqnAlgorithm', 'gn' ), Xv_km1_k, camera, yi );
%     
%         fprintf('residuo dh_dyi_newPar, %g\n',norm(J-Hii));

return



function a = dhrl_dy( Xv_km1_k, yi )

    rw = Xv_km1_k( 1:3 );
    Rrw = inv( q2r( Xv_km1_k( 4:7 ) ) );
    lambda = yi(6);
    phi = yi(5);
    theta = yi(4);
    
    dmi_dthetai = Rrw*[cos(phi)*cos(theta)  0   -cos(phi)*sin(theta)]';
    dmi_dphii = Rrw*[-sin(phi)*sin(theta)  -cos(phi)   -sin(phi)*cos(theta)]';
    
    a = [lambda*Rrw  dmi_dthetai dmi_dphii Rrw*(yi(1:3)-rw) ];
    
return



function Hi1 = dh_dxv( camera, Xv_km1_k, yi, zi )

    Hi1 = [ dh_drw( camera, Xv_km1_k, yi, zi )  dh_dqwr( camera, Xv_km1_k, yi, zi ) zeros( 2, 6 )];
    
%     % Verification, OK
%     [X,FVAL,EXITFLAG,OUTPUT,J] = fsolve( 'F_Test_dh_dxv_newPar', Xv_km1_k, ...
%             optimset( 'Display', 'off', ...
%             'NonlEqnAlgorithm', 'gn' ), yi, camera, Xv_km1_k);
% 
%     fprintf('residuo dh_dxv_newPar, %g\n',norm(J-Hi1));

return



function Hi12 = dh_dqwr( camera, Xv_km1_k, yi, zi )

    Hi12 = dh_dhrl( camera, Xv_km1_k, yi, zi ) * dhrl_dqwr( Xv_km1_k, yi );
    
return



function a = dhrl_dqwr( Xv_km1_k, yi )

    rw = Xv_km1_k( 1:3 );
    qwr = Xv_km1_k( 4:7 );
    lambda = yi(6);
    phi = yi(5);
    theta = yi(4);
    mi = [cos(phi)*sin(theta)   -sin(phi)  cos(phi)*cos(theta)]';
    
    a = dRq_times_a_by_dq( qconj(qwr), ((yi(1:3) - rw)*lambda + mi) )*dqbar_by_dq;
    
return



function Hi11 = dh_drw( camera, Xv_km1_k, yi, zi )

    Hi11 = dh_dhrl( camera, Xv_km1_k, yi, zi ) * dhrl_drw( Xv_km1_k, yi );

return



function a = dhrl_drw( Xv_km1_k, yi )
    
    % Verification, OK
    a = -( inv( q2r( Xv_km1_k(4:7) ) ) )*yi(6);

return



function a = dh_dhrl( camera, Xv_km1_k, yi, zi )

    a = dhd_dhu( camera, zi )*dhu_dhrl( camera, Xv_km1_k, yi );

return



function a = dhd_dhu( camera, zi_d )

    inv_a = jacob_undistor_fm( camera, zi_d );
    a=inv(inv_a);
    
%      % Verification, OK
%     [X,FVAL,EXITFLAG,OUTPUT,J] = fsolve( 'F_Test_dhd_dhu', zi_d, ...
%                                             optimset( 'Display', 'off', ...
%                                             'NonlEqnAlgorithm', 'gn' ), camera, zi_d );
%                                         
%     % fprintf('residuo dhd_dhu, %g\n',norm(J-a))
    
return
  
    
function a = dhu_dhrl( camera, Xv_km1_k, yi )
    
    f = camera.f;
    ku = 1/camera.dx;
    kv = 1/camera.dy;
    rw = Xv_km1_k( 1:3 );
    Rrw = inv(q2r( Xv_km1_k( 4:7 ) ));
    
    theta = yi(4);
    phi = yi(5);
    rho = yi(6);
    mi = [cos(phi)*sin(theta)   -sin(phi)  cos(phi)*cos(theta)]';
    
    hc = Rrw*( (yi(1:3) - rw)*rho + mi );
    hcx = hc(1);
    hcy = hc(2);
    hcz = hc(3);
    a = [+f*ku/(hcz)       0           -hcx*f*ku/(hcz^2);
        0               +f*kv/(hcz)    -hcy*f*kv/(hcz^2)];
    
%     % Verification, OK
%     [X,FVAL,EXITFLAG,OUTPUT,J] = fsolve( 'F_Test_dhu_dhrl', hc, ...
%                                             optimset( 'Display', 'off', ...
%                                             'NonlEqnAlgorithm', 'gn' ), camera, hc);
%                                         
%     fprintf('residuo dhd_dhrl, %g\n',norm(J-a));

return