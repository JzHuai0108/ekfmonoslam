function eul = dcm2euler(dcm)

pitch=asin(-dcm(3,1));  %pitch is assumed to be [-pi pi]. singular at pi. use ad-hoc methods to remedy this deficiency
roll=atan2(dcm(3,2), dcm(3,3));
heading = atan2(dcm(2,1), dcm(1,1));
eul=[roll;pitch;heading];