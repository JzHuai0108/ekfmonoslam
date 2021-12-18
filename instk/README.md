The instk package is downloaded from [here](http://www.instk.org/Tk/index.html).
Its documentation is at the [website](http://www.instk.org/Tk/Calibration/index.html).

Most files are identical to the original files from the instk website as of Dec 20 2021.
Only changes to three files are noteworthy.
In imu_err_defs_v000.m, extra IMU models are added.
In sys_metric_phipsi_v000.m, the IMU propagation model number is passed in, and Anav expression is corrected.
In AddIMUErr_v002.m, the IMU noise generation is updated.