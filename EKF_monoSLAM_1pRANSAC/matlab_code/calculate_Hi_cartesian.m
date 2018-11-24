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

function Hi = calculate_Hi_cartesian( Xv_km1_k, yi, camera, i, features_info )

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

% % Verification
% [X,FVAL,EXITFLAG,OUTPUT,J] = fsolve( 'F_Test_dh_dxv', Xv_km1_k, ...
%                                         optimset( 'Display', 'off', ...
%                                         'NonlEqnAlgorithm', 'gn' ), yi, camera, Xv_km1_k);
%                                         
% fprintf('residuo dh_dxv, %g\n',norm(J-Hi));

index_of_insertion = 13 + 3*sum(cartesian_features_index(1:i-1)) + 6*sum(inverse_depth_features_index(1:i-1))+1;
Hi(:,index_of_insertion:index_of_insertion+2) = dh_dy( camera, Xv_km1_k, yi, zi );

return



function Hii = dh_dy( camera, Xv_km1_k, yi, zi )

    Hii = dh_dhrl( camera, Xv_km1_k, yi, zi ) * dhrl_dy( Xv_km1_k );
    
%     % Verification
%     [X,FVAL,EXITFLAG,OUTPUT,J] = fsolve( 'F_Test_dh_dy', yi, ...
%                                             optimset( 'Display', 'off', ...
%                                             'NonlEqnAlgorithm', 'gn' ), camera, Xv_km1_k, yi);
%                                         
%     fprintf('residuo dh_dy, %g\n',norm(J-Hii));

return



function a = dhrl_dy( Xv_km1_k )

    a = inv( q2r( Xv_km1_k( 4:7 ) ) );
    
return



function Hi1 = dh_dxv( camera, Xv_km1_k, yi, zi )

    Hi1 = [ dh_drw( camera, Xv_km1_k, yi, zi )  dh_dqwr( camera, Xv_km1_k, yi, zi ) zeros( 2, 6 )];

return



function Hi12 = dh_dqwr( camera, Xv_km1_k, yi, zi )

    Hi12 = dh_dhrl( camera, Xv_km1_k, yi, zi ) * dhrl_dqwr( Xv_km1_k, yi );
    
return



function a = dhrl_dqwr( Xv_km1_k, yi )

    a = dRq_times_a_by_dq( qconj(Xv_km1_k( 4:7 )), (yi - Xv_km1_k( 1:3 )) )*dqbar_by_dq;
    
return



function Hi11 = dh_drw( camera, Xv_km1_k, yi, zi )

    Hi11 = dh_dhrl( camera, Xv_km1_k, yi, zi ) * dhrl_drw( Xv_km1_k );

return



function a = dhrl_drw( Xv_km1_k )

    a = -( inv( q2r( Xv_km1_k(4:7) ) ) );

return



function a = dh_dhrl( camera, Xv_km1_k, yi, zi )

    a = dhd_dhu( camera, zi )*dhu_dhrl( camera, Xv_km1_k, yi );

return



function a = dhd_dhu( camera, zi_d )

    inv_a = jacob_undistor_fm( camera, zi_d );
    a=inv(inv_a);
    
%      % Verification, OK
%     [X,FVAL,EXITFLAG,OUTPUT,J] = fsolve( 'F_Test_dhd_dhu', zi_u, ...
%                                             optimset( 'Display', 'off', ...
%                                             'NonlEqnAlgorithm', 'gn' ), camera, zi_u );
%                                         
%     fprintf('residuo dhd_dhu, %g\n',norm(J-a))
    
return
  
    
function a = dhu_dhrl( camera, Xv_km1_k, yi )
    
    f = camera.f;
    ku = 1/camera.dx;
    kv = 1/camera.dy;
    rw = Xv_km1_k( 1:3 );
    Rrw = inv(q2r( Xv_km1_k( 4:7 ) ));
    hrl = Rrw*( yi - rw );
    hrlx = hrl(1);
    hrly = hrl(2);
    hrlz = hrl(3);
    a = [f*ku/(hrlz)       0           -hrlx*f*ku/(hrlz^2);
        0               f*kv/(hrlz)    -hrly*f*kv/(hrlz^2)];
    
% %     Verification
%     [X,FVAL,EXITFLAG,OUTPUT,J] = fsolve( 'F_Test_dhu_dhrl', hrl, ...
%                                             optimset( 'Display', 'off', ...
%                                             'NonlEqnAlgorithm', 'gn' ), camera, hrl);
%                                         
%     fprintf('residuo dhd_dhrl, %g\n',norm(J-a));

return