Javier Civera, jcivera@unizar.es, 25/june/2012

Quick-and-dirty code for converting the camera calibration model of the popular Bouguet's Calibration Toolbox [1] to the camera calibration model used in our EKF Monocular SLAM Matlab code [2]. 

The code basically does a non-linear optimization, minimizing the error between some image points aplying Bouguet's calibration model and the one we use. 

In order to make it work, just copy your calibration results (the .mat file) from Bouguet's toolbox in the directory, go to the .m file Bouguet2Tsai.m, change the name of the calibration results file in line 5, and run it. The result should be the calibration parameters you should use in the EKF monocular SLAM code. You have details of the output in the Bouguet2Tsai.m file. Also, if you run the code, you can see an example of the conversion done.

[1] http://www.vision.caltech.edu/bouguetj/calib_doc/

[2] http://webdiis.unizar.es/~jcivera/code/1p-ransac-ekf-monoslam.html
