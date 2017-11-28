# ekfmonoslam

SLAM based on EKF using a monocular camera, optionally an IMU, and GPS data. This package has been tested in matlab 2012 on 64bit Windows 8.1. 

Before running sample tests, a few packages need to be downloaded. To begin with, create a folder to contain all the related packages, e.g., matlab_ws.

# 1. Dependencies

## 1.1 1-Point RANSAC Inverse Depth EKF Monocular SLAM v1.01

The EKF of ekfmonoslam is based on that in 1-Point RANSAC Inverse Depth EKF Monocular SLAM v1.01. You can download the package from http://webdiis.unizar.es/~jcivera/code/1p-ransac-ekf-monoslam.html. Extract the package, and put the EKF_monoSLAM_1pRANSAC folder into matlab_ws. Note this package contains an example sequence. Change the name of the subfolder /EKF_monoSLAM_1pRANSAC/matlab_code/@ekf_filter into /EKF_monoSLAM_1pRANSAC/matlab_code/ekf_filter

## 1.2 Camera calibration toolbox for matlab

Ekfmonoslam used a couple of pinhole camera projection functions from camera calibration toolbox. From this web page http://www.vision.caltech.edu/bouguetj/calib_doc/, search for "download page", then download the package and extract and put it into matlab_ws, too.

## 1.3 MexOpenCV and OpenCV

Ekfmonoslam used Good feature to track and Lucas-Kanade tracker in OpenCV. MexOpenCV can be forked or downloaded from https://github.com/kyamagu/mexopencv. Also, extract and put mexopencv into matlab_ws folder. Compile the mex files according to the instructions on https://github.com/kyamagu/mexopencv.

## 1.4 INS toolkit (included)

ekfmonoslam also depends on the INS toolkit which can be found at http://www.instk.org/. Since I made much changes to instk, I included it in ekfmonoslam package.

## 1.5 VoiceBox

Mike Brooke's VoiceBox are used for rotation transformations. It can be downloaded from http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html. Again, put it in matlab_ws folder.

# 2. Build

Download the ekfmonoslam package, put all folders within ekfmonoslam into matlab_ws. Because there is a EKF_monoSLAM_1pRANSAC folder in this package which contains patches to 1-Point RANSAC Inverse Depth EKF Monocular SLAM v1.01, you may be prompted to replace or merge files, replace these files already in matlab_ws.

# 3. Test with data

A test data can be retrieved from [here](https://pan.baidu.com/s/1c1IdiQO). Once downloaded, put the files and a folder within it to the matlab_ws/data folder. To run a test, in matlab navigate to the folder of matlab_ws/ekfmonocularslamv02, then run mainslamhuai_v002.m. Its first experiment should execute by default.

# 4. Simulate noisy IMU data given ground truth poses

Suppose the ground truth SE(3) poses (3D position and 3DOF attitude) are given at a regular interval, e.g. 0.1s, the noisy IMU data at a particular frequency, e.g., 100Hz can be generated in two steps: 1, generate the continuous trajectory by interpolating the existing poses, differentiating the continuous trajectory will give angular rates and linear accelerations, for details, see [1]; 2, add noise to these true IMU samples, given error specs of a particular IMU model.

In practice, first build and run the program in the folder SE3Interpolation. Two test cases are provided in SE3Interpolation, for both KITTI seq 00 and Tsukuba CG stereo dataset. To run it, you need to download at least one of these two datasets. This step will produce the true IMU samples. Then to add noise, call testAddImuError.m in utilities/tests folder.


# Reference

[1] Jianzhu Huai, Charles K. Toth and Dorota A. Grejner-Brzezinska: Stereo-inertial odometry using nonlinear optimization. Proceedings of the 28th International Technical Meeting of The Satellite Division of the Institute of Navigation (ION GNSS+ 2015) 2015.
