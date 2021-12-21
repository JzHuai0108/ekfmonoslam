# ekfmonoslam

SLAM based on EKF using a monocular camera, optionally an IMU, and GPS data. This package has been tested in matlab 2012 on 64bit Windows 8.1. 

# 1. Dependencies

## 1.1 1-Point RANSAC Inverse Depth EKF Monocular SLAM v1.01 (included)

The EKF of ekfmonoslam is based on that in 1-Point RANSAC Inverse Depth [EKF Monocular SLAM v1.01](http://webdiis.unizar.es/~jcivera/code/1p-ransac-ekf-monoslam.html). A few minor changes have been made to several files of its scripts, so the package is included to simplify the installation procedure.

## 1.2 Camera calibration toolbox for matlab (included)

Ekfmonoslam used a couple of pinhole camera projection functions from the [camera calibration toolbox package](http://www.vision.caltech.edu/bouguetj/calib_doc/download/toolbox_calib.zip).

## 1.3 mexopencv (included)

Ekfmonoslam used Good feature to track and Lucas-Kanade tracker in [mexopencv](https://github.com/kyamagu/mexopencv/). Mex files built from either mexopencv v3.4.1 or v2.4 have been shown to support ekfmonoslam. For convenience, necessary compiled mex files are included for [mexopencv v2.4](https://github.com/kyamagu/mexopencv/tree/v2.4) built in maplab 2012b on Windows 7 with OpenCV 2.4.6 x64 vc10 release libs. These mex files worked normally under maplab 2012b and 2018b. If mexopencv is to be built anew, please refer to the instuctions [here](https://github.com/kyamagu/mexopencv/).

## 1.4 INS toolkit (included)

ekfmonoslam also depends on the INS toolkit which can be found at http://www.instk.org/. Since I made much changes to instk, I included it in ekfmonoslam package.

## 1.5 VoiceBox (included)

Mike Brooke's [VoiceBox](http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html) is used for rotation transformations.


# 2. Build

Download the ekfmonoslam package.

# 3. Test with data

A test data, "20130808", can be retrieved from [here](https://pan.baidu.com/s/1c1IdiQO). Once downloaded, put it to the ekfmonoslam/data folder.
(1) To test GPS/IMU integration, open the folder ekfmonoslam/tests in matlab, then run Test_EKF_filters.m.

(2) To test GPS/IMU/Monocular fusion, run mainslam.m.

(3) Unscented Kalman filter for GPS/IMU/Monocular camera fusion, 
refer to gssm_cam_ins_gps.m and cam_ins_gps.m. 
But this was very experimental.

# 4. Simulate noisy IMU data given ground truth poses

Suppose the ground truth SE(3) poses (3D position and 3DOF attitude) are given at a regular interval, e.g. 0.1s, the noisy IMU data at a particular frequency, e.g., 100Hz can be generated in two steps: 1, generate the continuous trajectory by interpolating the existing poses, differentiating the continuous trajectory will give angular rates and linear accelerations, for details, see [1]; 2, add noise to these true IMU samples, given error specs of a particular IMU model.

In practice, first build and run the program in the folder SE3Interpolation. Two test cases are provided in SE3Interpolation, for both KITTI seq 00 and Tsukuba CG stereo dataset. To run it, you need to download at least one of these two datasets. This step will produce the true IMU samples. Then to add noise, call testAddImuError.m in tests folder.


# Reference

[1] Jianzhu Huai, Charles K. Toth and Dorota A. Grejner-Brzezinska: Stereo-inertial odometry using nonlinear optimization. Proceedings of the 28th International Technical Meeting of The Satellite Division of the Institute of Navigation (ION GNSS+ 2015) 2015.
