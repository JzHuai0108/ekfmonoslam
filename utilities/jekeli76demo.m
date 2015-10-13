% demonstrate the kalman filter in Jekeli, 2000, ch7.6 example p.227
experim=1; % 0 for position measurement, 1 for velocity measurements
Ndata=100;
initX=0;
initV=0;
deltat=1;
atrue=2;
biastrue=1;
Qpk=0.09;
Qvk=0.01;
procnoiseq=0.01;
time=(1:Ndata)*deltat;
Xtrue=zeros(2,Ndata);
Xtrue(1,:)=initX+initV*time+0.5*atrue*time.^2;
Xtrue(2,:)=initV+atrue*time;
Ypos=Xtrue(1,:)+randn(1, Ndata)*chol(Qpk);% measurements of position from t1 to tN
Yvel=Xtrue(2,:)+randn(1, Ndata)*chol(Qvk);% measurements of velocity from t1 to tN
deltav=(atrue+biastrue)*deltat+randn(1, Ndata)*chol(procnoiseq*deltat);
aindicate=deltav/deltat;% odometry input, from t1 to tN
filename='deleteme.txt';
filter = EKF_filter_jekeli76(initX, initV, deltat,procnoiseq );
ffilres=fopen(filename,'wt');
fprintf(ffilres, '%% time, estimated, vel, pos, bias, vel cov, pos cov, bias cov\n');
filter.SaveToFile(0, ffilres);

for rider=1:Ndata
    % predict states and covariance
    filter.ffun_state(aindicate(rider));
    filter.ffun_covariance();
    % grab a observation and correct the states
    if(mod(rider,3)==0)
       
        switch(experim)
            case 0
                 predict=filter.xi;
                measure=Ypos(rider);
                      H=[0, 1,0];
        R=Qpk;
            case 1
                 predict=filter.vi;
                measure=Yvel(rider);
                      H=[1, 0,0];
        R=Qvk;
        end
  
        filter.correctstates(predict,measure, H,R);
    end
    filter.SaveToFile(rider*deltat, ffilres);
end
res=load(filename);
figure(1)
plot(res(:, 1), res(:, 4), '-r+');
hold on
plot(res(:, 1), res(:, 5),'-kd');
hold on
plot(res(:, 1), res(:, 6),'-gs');
hold on
plot(res(:, 1), res(:, 7),'-mo');
hold off
legend('bias', '\sigma_v', '\sigma_x', '\sigma_b');
figure (2)
plot(res(:,1), res(:,2)-[0, Xtrue(2,:)]','-ks');
hold on
plot(res(:,1), res(:,3)-[0, Xtrue(1,:)]','-r^');
legend('velocity','position')