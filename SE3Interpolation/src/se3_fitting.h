#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "sophus/se3.hpp"

// interpolate IMU data given poses and their uniform timestamps
// input: q02n, times: q_0^w, q_1^w, ..., q_n^w; N=n+1 poses and their timestamps; also outputFreq
// output: samplePoses: sampled poses; samples: output timestamps, acceleration of sensor by combined force in world frame,
// and angular rate of sensor w.r.t world frame represented in sensor frame, and velocity of sensor in world frame
template<class Scalar>
void InterpolateIMUData(const std::vector<Sophus::SE3Group<Scalar> > &q02n,const std::vector<Scalar>& times, const Scalar outputFreq,
                        std::vector<Eigen::Matrix<Scalar, 4,4 > >& samplePoses, std::vector<Eigen::Matrix<Scalar, 10, 1> >& samples)
{
    typedef Sophus::SO3Group<Scalar> SO3Type;
    typedef Sophus::SE3Group<Scalar> SE3Type;
    typedef typename Sophus::SE3Group<Scalar>::Tangent Tangent;

    std::cout<<"Assigning control points"<<std::endl;
    int lineNum=q02n.size();
    std::vector<SE3Type> bm32nm1(lineNum+2); // b_-3^w, b_-2^w, ..., b_(n-1)^w
    std::vector<Tangent> Omegam22nm1(lineNum+1); // $\Omega_-2, \Omega_-1, ..., \Omega_n-1$ where $\Omega_j=log((b_j-1)\b_j)$
// assume initial velocity is zero
//    Omegam22nm1[1]=SE3Type::log(q02n[0].inverse()*q02n[1]); //\Omega_-1
//    Omegam22nm1[0]=SE3Type::vee(-SE3Type::exp(Omegam22nm1[1]/6).matrix()*SE3Type::hat(Omegam22nm1[1])*
//            SE3Type::exp(-Omegam22nm1[1]/6).matrix());
//    bm32nm1[1]=q02n[0];
//    bm32nm1[2]=q02n[1];
//    bm32nm1[0]=bm32nm1[1]*SE3Type::exp(-Omegam22nm1[0]);
// or assume first three poses have identical difference
    bm32nm1[1]=q02n[0];
    bm32nm1[2]=q02n[1];
    bm32nm1[0]=bm32nm1[1]*bm32nm1[2].inverse()*bm32nm1[1];
    Omegam22nm1[0]=SE3Type::log(bm32nm1[0].inverse()*bm32nm1[1]);
    Omegam22nm1[1]=SE3Type::log(q02n[0].inverse()*q02n[1]); //\Omega_-1

    for (int i=3; i<lineNum+1; ++i)
    {
        bm32nm1[i]=q02n[i-1];
        Omegam22nm1[i-1]=SE3Type::log(bm32nm1[i-1].inverse()*bm32nm1[i]);
    }
    bm32nm1[lineNum+1]=q02n[lineNum-1]*q02n[lineNum-2].inverse()*q02n[lineNum-1];
    Omegam22nm1[lineNum]=SE3Type::log(bm32nm1[lineNum].inverse()*bm32nm1[lineNum+1]);

    std::cout<<"take derivatives to compute acceleration and angular rate"<<std::endl;
    int dataCount=floor((*(times.rbegin())-1e-6-times[0])*outputFreq)+1; // how many output data, from t_0 up to close to t_n
    samplePoses.resize(dataCount);
    samples.resize(dataCount); // output timestamps, acceleration of sensor in world frame,
    // and angular rate of sensor w.r.t world frame represented in sensor frame

    Eigen::Matrix<Scalar, 4,4> sixC; // six times C matrix
    sixC<<6, 0, 0, 0,
            5, 3, -3, 1,
            1, 3, 3, -2,
            0, 0, 0, 1;
    Scalar timestamp, Deltat, ut;
    Eigen::Matrix<Scalar, 4,1> utprod, tildeBs, dotTildeBs, ddotTildeBs;
    std::vector<SE3Type> tripleA(3);//A_1, A_2, A_3
    std::vector<Eigen::Matrix<Scalar, 4,4> > dotDdotAs(6);
    //$\dot{A_1}, \dot{A_2}, \dot{A_3}, \ddot{A_1}, \ddot{A_2}, \ddot{A_3}$
    // where $p(t)=b_{i-3}*A_1*A_2*A_3$ for $t\in[t_i, t_{i+1})$
    SE3Type Ts2w; //T_s^w
    std::vector<Eigen::Matrix<Scalar, 4,4> > dotDdotTs(2); //$\dot{T_s^w}, \ddot{T_s^w}$
    int tickIndex=0; // where is a timestamp in times, s.t. $timestamp\in[t_{tickIndex}, t_{tickIndex+1})$
    for (int i=0; i<dataCount;++i){
        timestamp=times[0]+i/outputFreq;
        samples[i][0]=timestamp;
        if(timestamp>=times[tickIndex+1])
            tickIndex=tickIndex+1;
        assert(timestamp<times[tickIndex+1]);

        Deltat=times[tickIndex+1]-times[tickIndex];
        ut=(timestamp-times[tickIndex])/Deltat;
        utprod<<1, ut, ut*ut, ut*ut*ut;
        tildeBs=sixC*utprod/6;
        utprod<<0, 1, 2*ut, 3*ut*ut;
        dotTildeBs=sixC*utprod/(6*Deltat);
        utprod<<0, 0, 2, 6*ut;
        ddotTildeBs=sixC*utprod/(6*Deltat*Deltat);

        tripleA[0]=SE3Type::exp(Omegam22nm1[tickIndex]*tildeBs[1]);
        tripleA[1]=SE3Type::exp(Omegam22nm1[tickIndex+1]*tildeBs[2]);
        tripleA[2]=SE3Type::exp(Omegam22nm1[tickIndex+2]*tildeBs[3]);
        dotDdotAs[0]=tripleA[0].matrix()*SE3Type::hat(Omegam22nm1[tickIndex])*dotTildeBs[1];
        dotDdotAs[1]=tripleA[1].matrix()*SE3Type::hat(Omegam22nm1[tickIndex+1])*dotTildeBs[2];
        dotDdotAs[2]=tripleA[2].matrix()*SE3Type::hat(Omegam22nm1[tickIndex+2])*dotTildeBs[3];
        dotDdotAs[3]=tripleA[0].matrix()*SE3Type::hat(Omegam22nm1[tickIndex])*ddotTildeBs[1]+
                dotDdotAs[0]*SE3Type::hat(Omegam22nm1[tickIndex])*dotTildeBs[1];
        dotDdotAs[4]=tripleA[1].matrix()*SE3Type::hat(Omegam22nm1[tickIndex+1])*ddotTildeBs[2]+
                dotDdotAs[1]*SE3Type::hat(Omegam22nm1[tickIndex+1])*dotTildeBs[2];
        dotDdotAs[5]=tripleA[2].matrix()*SE3Type::hat(Omegam22nm1[tickIndex+2])*ddotTildeBs[3]+
                dotDdotAs[2]*SE3Type::hat(Omegam22nm1[tickIndex+2])*dotTildeBs[3];

        Ts2w=bm32nm1[tickIndex]*tripleA[0]*tripleA[1]*tripleA[2];
        dotDdotTs[0]=bm32nm1[tickIndex].matrix()*dotDdotAs[0]*(tripleA[1]*tripleA[2]).matrix()+
                (bm32nm1[tickIndex]*tripleA[0]).matrix()*dotDdotAs[1]*tripleA[2].matrix()+
                (bm32nm1[tickIndex]*tripleA[0]*tripleA[1]).matrix()*dotDdotAs[2];

        dotDdotTs[1]=bm32nm1[tickIndex].matrix()*dotDdotAs[3]*(tripleA[1]*tripleA[2]).matrix()+
                (bm32nm1[tickIndex]*tripleA[0]).matrix()*dotDdotAs[4]*tripleA[2].matrix()+
                (bm32nm1[tickIndex]*tripleA[0]*tripleA[1]).matrix()*dotDdotAs[5]+
                2*bm32nm1[tickIndex].matrix()*(dotDdotAs[0]*dotDdotAs[1]*tripleA[2].matrix()+
                tripleA[0].matrix()*dotDdotAs[1]*dotDdotAs[2]+dotDdotAs[0]*tripleA[1].matrix()*dotDdotAs[2]);

        samplePoses[i]=Ts2w.matrix();
        samples[i].segment(1,3)=dotDdotTs[1].col(3).head(3);//$a_s^w$        
        samples[i].segment(4,3)=SO3Type::vee(Ts2w.rotationMatrix().transpose()*dotDdotTs[0].topLeftCorner(3,3));//$\omega_{ws}^s$
        samples[i].tail(3)=dotDdotTs[0].col(3).head(3);//$v_s^w$        
    }
}

