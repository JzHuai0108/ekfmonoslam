s = RandStream('mt19937ar','Seed',1);
reset(s);
acc=randn(s,3,25);
gyro=randn(s,3,25);

[vinc, vcorr]=sculling_v000(gyro, acc, 4, 5);
[vinc1, vcorr1, ainc1, acorr1]=sculling_coning_v000(gyro, acc, 4, 1);

vcorr-vcorr1
