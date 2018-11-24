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

function dqomegadt_by_domegaRES=dqomegadt_by_domega(omega, delta_t)

%  // Calculate commonly used Jacobian part dq(omega * delta_t) by domega

  %// Modulus
  omegamod = norm(omega);

  %// Use generic ancillary functions to calculate components of Jacobian
  dqomegadt_by_domegaRES(1, 1) = dq0_by_domegaA(omega(1), omegamod, delta_t);
  dqomegadt_by_domegaRES(1, 2) = dq0_by_domegaA(omega(2), omegamod, delta_t);
  dqomegadt_by_domegaRES(1, 3) = dq0_by_domegaA(omega(3), omegamod, delta_t);
  dqomegadt_by_domegaRES(2, 1) = dqA_by_domegaA(omega(1), omegamod, delta_t);
  dqomegadt_by_domegaRES(2, 2) = dqA_by_domegaB(omega(1), omega(2), omegamod, delta_t);
  dqomegadt_by_domegaRES(2, 3) = dqA_by_domegaB(omega(1), omega(3), omegamod, delta_t);
  dqomegadt_by_domegaRES(3, 1) = dqA_by_domegaB(omega(2), omega(1), omegamod, delta_t);
  dqomegadt_by_domegaRES(3, 2) = dqA_by_domegaA(omega(2), omegamod, delta_t);
  dqomegadt_by_domegaRES(3, 3) = dqA_by_domegaB(omega(2), omega(3), omegamod, delta_t);
  dqomegadt_by_domegaRES(4, 1) = dqA_by_domegaB(omega(3), omega(1), omegamod, delta_t);
  dqomegadt_by_domegaRES(4, 2) = dqA_by_domegaB(omega(3), omega(2), omegamod, delta_t);
  dqomegadt_by_domegaRES(4, 3) = dqA_by_domegaA(omega(3), omegamod, delta_t);
  
 return


% // Ancillary functions: calculate parts of Jacobian dq_by_domega
% // which are repeatable due to symmetry.
% // Here omegaA is one of omegax, omegay, omegaz
% // omegaB, omegaC are the other two
% // And similarly with qA, qB, qC

function dq0_by_domegaARES=dq0_by_domegaA(omegaA, omega, delta_t)

  dq0_by_domegaARES=(-delta_t / 2.0) * (omegaA / omega) * sin(omega * delta_t / 2.0);
  

function dqA_by_domegaARES=dqA_by_domegaA(omegaA, omega, delta_t)
  dqA_by_domegaARES=(delta_t / 2.0) * omegaA * omegaA / (omega * omega) ...
                    * cos(omega * delta_t / 2.0) ...
                    + (1.0 / omega) * (1.0 - omegaA * omegaA / (omega * omega))...
                    * sin(omega * delta_t / 2.0);

return

function dqA_by_domegaBRES=dqA_by_domegaB(omegaA, omegaB, omega, delta_t)

  dqA_by_domegaBRES=(omegaA * omegaB / (omega * omega)) * ...
                    ( (delta_t / 2.0) * cos(omega * delta_t / 2.0) ...
                    - (1.0 / omega) * sin(omega * delta_t / 2.0) );
return    
            
