=============================================================
1-Point RANSAC Inverse Depth EKF-Based Structure from Motion. 
Version: 1.0.
=============================================================

Copyright (C) 2010 Javier Civera and J. M. M. Montiel
Universidad de Zaragoza, Zaragoza, Spain.
 
Matlab code for EKF-Based Structure from Motion / EKF SLAM from a monocular sequence. 

------------------------------------------------------------
License
------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation. Read http://www.gnu.org/copyleft/gpl.html (included here in GPL.txt) for details. In short, that means:

1) You can freely use and modify this code.

2) If you want to distribute code based on this one, it has to be done under the GNU GPL License.

If you use this code for academic work, please reference:

[1]   Javier Civera, Oscar G. Grasa, Andrew J. Davison, J. M. M. Montiel,
      1-Point RANSAC for EKF Filtering: Application to Real-Time Structure from Motion and Visual Odometry,
      to appear in Journal of Field Robotics, October 2010.

For commercial uses, please contact us.

------------------------------------------------------------
Quick start
------------------------------------------------------------

After unzipping:

1) Go to matlab_code directory
2) Run mono_slam.m

The zipped file includes an example sequence, so the code should run as it is.
The code has been tested from Matlab 6.5 to Matlab 7.9. 

It should not depend on any specific Matlab toolbox.

The directory matlab_code/figures contains a video displaying what the code should show.

------------------------------------------------------------
Support
------------------------------------------------------------

Questions, comments, suggestions, discussions and bug reports are welcomed. Please, mail to jcivera@unizar.es

------------------------------------------------------------
Dependencies
------------------------------------------------------------

In this matlab code we are using the FAST corner detector by Edward Rosten and Tom Drummond [2], downloaded from http://mi.eng.cam.ac.uk/~er258/work/fast.html.

[2]   Edward Rosten and Tom Drummond,
      Machine learning for high-speed corner detection,
      European Conference on Computer Vision, May 2006.

We are also using some functions from the Robotics Toolbox for Matlab by Peter Corke [3], downloaded from http://www.petercorke.com/Robotics_Toolbox.html.

[3]   P.I. Corke
      MATLAB toolboxes: robotics and vision for students and teachers
      IEEE Robotics and Automation Magazine, Volume 14(4), December 2007, pp. 16-17

=============================================================
Versions history
=============================================================

1.0
---
- Initial release
- Known things to fix in later versions:
  a) Only Tsai distortion model available
  b) Robocentric not implemented yet.
  c) EKF filter class rather messy.

This code has been developed for research purposes. Please, do expect to find some illustrative examples of bad programming practices :)
