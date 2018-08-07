#ifndef SE3_FITTING_H_
#define SE3_FITTING_H_

#include <iostream>
#include <vector>

#include <Eigen/Core>

#include "sophus/se3.hpp"

/**
 *@brief interpolate IMU data given control poses and their uniform timestamps
 * @param q02n, nominal trajecotry poses,i.e., control points, q_0^w, q_1^w, ..., q_n^w; 
 * N=n+1 sensor to world transforms
 * for interpolation, one pose is added at both ends of the array of q02n, making its size n+3.
 * The two poses are added assuming constant velocity at the start and the end.
 * @param times, N timestamps, assume equally spaced in time
 * @param outputFreq, output frequency of true inertial data
 * @param samplePoses, output sampled sensor to world transforms
 * @param samples, output inertial samples. Each entry: timestamps, acceleration of sensor
 * by the combined force (specific forces and gravity) in the world frame if gw is not specified,
 * otherwise, acceleration of sensor by the specific force in the sensor frame,
 * and angular rate of sensor w.r.t world frame represented in the sensor frame, 
 * and velocity of sensor in world frame
 * output timestamps are times[0], times[0]+1/outputFreq, times[0]+2/outputFreq, ...
*/
template<class Scalar>
void InterpolateIMUData(const std::vector<Sophus::SE3Group<Scalar>, Eigen::aligned_allocator<Sophus::SE3Group<Scalar> > > &q02n,
                        const std::vector<Scalar>& times,
                        const Scalar outputFreq,
                        std::vector<Eigen::Matrix<Scalar, 4, 4>, Eigen::aligned_allocator<Eigen::Matrix<Scalar, 4, 4> > >& samplePoses,
                        std::vector<Eigen::Matrix<Scalar, 10, 1>, Eigen::aligned_allocator<Eigen::Matrix<Scalar, 10, 1> > >& samples, 
                        const Eigen::Matrix<Scalar, 3, 1> gw = Eigen::Matrix<Scalar, 3, 1>::Zero())
{
    bool outputaccelinw = true;
    if (!gw.norm() == 0)
      outputaccelinw = false;
    typedef Sophus::SO3Group<Scalar> SO3Type;
    typedef Sophus::SE3Group<Scalar> SE3Type;
    typedef typename Sophus::SE3Group<Scalar>::Tangent Tangent;

    std::cout<<"Assigning control points"<<std::endl;
    int lineNum=q02n.size();
    std::vector<SE3Type, Eigen::aligned_allocator<SE3Type> > bm12np1(lineNum+2); // b_-1^w, b_0^w, ..., b_(n+1)^w
    std::vector<Tangent, Eigen::aligned_allocator<Tangent> > Omega02np1(lineNum+1); // $\Omega_0, \Omega_1, ..., \Omega_n+1$ where $\Omega_j=log((b_j-1)\b_j)$
// assume initial velocity is zero, often cause jumps in acceleration, not recommended
//    Omega02np1[1]=SE3Type::log(q02n[0].inverse()*q02n[1]); //\Omega_1
//    Omega02np1[0]=SE3Type::vee(-SE3Type::exp(Omega02np1[1]/6).matrix()*SE3Type::hat(Omega02np1[1])*
//            SE3Type::exp(-Omega02np1[1]/6).matrix());
//    bm12np1[1]=q02n[0];
//    bm12np1[2]=q02n[1];
//    bm12np1[0]=bm12np1[1]*SE3Type::exp(-Omega02np1[0]);
// or assume first three poses have identical difference
    bm12np1[1]=q02n[0];
    bm12np1[2]=q02n[1];
    bm12np1[0]=bm12np1[1]*bm12np1[2].inverse()*bm12np1[1];
    Omega02np1[0]=SE3Type::log(bm12np1[0].inverse()*bm12np1[1]);
    Omega02np1[1]=SE3Type::log(q02n[0].inverse()*q02n[1]); //\Omega_1

    for (int i=3; i<lineNum+1; ++i)
    {
        bm12np1[i]=q02n[i-1];
        Omega02np1[i-1]=SE3Type::log(bm12np1[i-1].inverse()*bm12np1[i]);
    }
    bm12np1[lineNum+1]=q02n[lineNum-1]*q02n[lineNum-2].inverse()*q02n[lineNum-1];
    Omega02np1[lineNum]=SE3Type::log(bm12np1[lineNum].inverse()*bm12np1[lineNum+1]);

    std::cout<<"Take derivatives to compute acceleration and angular rate"<<std::endl;
    int dataCount=floor((*(times.rbegin())-1e-6-times[0])*outputFreq)+1; // how many output data, from t_0 up to close to t_n
    samplePoses.resize(dataCount);
    samples.resize(dataCount);
    Eigen::Matrix<Scalar, 4,4> sixC; // six times C matrix
    sixC<<6, 0, 0, 0,
            5, 3, -3, 1,
            1, 3, 3, -2,
            0, 0, 0, 1;
    Scalar timestamp, Deltat, ut;
    Eigen::Matrix<Scalar, 4,1> utprod, tildeBs, dotTildeBs, ddotTildeBs;
    std::vector<SE3Type, Eigen::aligned_allocator<SE3Type> > tripleA(3);//A_1, A_2, A_3
    std::vector<Eigen::Matrix<Scalar, 4, 4>, Eigen::aligned_allocator<Eigen::Matrix<Scalar, 4, 4> > > dotDdotAs(6);
    //$\dot{A_1}, \dot{A_2}, \dot{A_3}, \ddot{A_1}, \ddot{A_2}, \ddot{A_3}$
    // where $p(t)=b_{i-3}*A_1*A_2*A_3$ for $t\in[t_i, t_{i+1})$
    SE3Type Ts2w; //T_s^w
    std::vector<Eigen::Matrix<Scalar, 4, 4>,
            Eigen::aligned_allocator<Eigen::Matrix<Scalar, 4, 4>> > dotDdotTs(2); //$\dot{T_s^w}, \ddot{T_s^w}$
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

        tripleA[0]=SE3Type::exp(Omega02np1[tickIndex]*tildeBs[1]);
        tripleA[1]=SE3Type::exp(Omega02np1[tickIndex+1]*tildeBs[2]);
        tripleA[2]=SE3Type::exp(Omega02np1[tickIndex+2]*tildeBs[3]);
        dotDdotAs[0]=tripleA[0].matrix()*SE3Type::hat(Omega02np1[tickIndex])*dotTildeBs[1];
        dotDdotAs[1]=tripleA[1].matrix()*SE3Type::hat(Omega02np1[tickIndex+1])*dotTildeBs[2];
        dotDdotAs[2]=tripleA[2].matrix()*SE3Type::hat(Omega02np1[tickIndex+2])*dotTildeBs[3];
        dotDdotAs[3]=tripleA[0].matrix()*SE3Type::hat(Omega02np1[tickIndex])*ddotTildeBs[1]+
                dotDdotAs[0]*SE3Type::hat(Omega02np1[tickIndex])*dotTildeBs[1];
        dotDdotAs[4]=tripleA[1].matrix()*SE3Type::hat(Omega02np1[tickIndex+1])*ddotTildeBs[2]+
                dotDdotAs[1]*SE3Type::hat(Omega02np1[tickIndex+1])*dotTildeBs[2];
        dotDdotAs[5]=tripleA[2].matrix()*SE3Type::hat(Omega02np1[tickIndex+2])*ddotTildeBs[3]+
                dotDdotAs[2]*SE3Type::hat(Omega02np1[tickIndex+2])*dotTildeBs[3];

        Ts2w=bm12np1[tickIndex]*tripleA[0]*tripleA[1]*tripleA[2];
        dotDdotTs[0]=bm12np1[tickIndex].matrix()*dotDdotAs[0]*(tripleA[1]*tripleA[2]).matrix()+
                (bm12np1[tickIndex]*tripleA[0]).matrix()*dotDdotAs[1]*tripleA[2].matrix()+
                (bm12np1[tickIndex]*tripleA[0]*tripleA[1]).matrix()*dotDdotAs[2];

        dotDdotTs[1]=bm12np1[tickIndex].matrix()*dotDdotAs[3]*(tripleA[1]*tripleA[2]).matrix()+
                (bm12np1[tickIndex]*tripleA[0]).matrix()*dotDdotAs[4]*tripleA[2].matrix()+
                (bm12np1[tickIndex]*tripleA[0]*tripleA[1]).matrix()*dotDdotAs[5]+
                2*bm12np1[tickIndex].matrix()*(dotDdotAs[0]*dotDdotAs[1]*tripleA[2].matrix()+
                tripleA[0].matrix()*dotDdotAs[1]*dotDdotAs[2]+dotDdotAs[0]*tripleA[1].matrix()*dotDdotAs[2]);

        samplePoses[i]=Ts2w.matrix();
        if (outputaccelinw)
            samples[i].segment(1,3)=dotDdotTs[1].col(3).head(3);//$a_s^w$
        else
            samples[i].segment(1,3)=Ts2w.unit_quaternion().inverse()*(dotDdotTs[1].col(3).head(3) - gw);//$a_m^s$
        samples[i].segment(4,3)=SO3Type::vee(Ts2w.rotationMatrix().transpose()*dotDdotTs[0].topLeftCorner(3,3));//$\omega_{ws}^s$
        samples[i].tail(3)=dotDdotTs[0].col(3).head(3);//$v_s^w$
    }
}

//input: lat and long, height is not needed
//output: Ce2n
inline Eigen::Matrix3d llh2dcm( Eigen::Vector3d &llh)
{
    double sL = sin(llh[0]);
    double cL = cos(llh[0]);
    double sl = sin(llh[1]);
    double cl = cos(llh[1]);

    Eigen::Matrix3d Ce2n;
    Ce2n<< -sL * cl, -sL * sl, cL ,  -sl, cl, 0 , -cL * cl, -cL * sl, -sL;
    return Ce2n;
}

//output rotation matrix about 3 axis for rad in radians
// if we rotate b frame about its 3 axis by rad into a frame, then =x^a*R3(rad)*x^b
inline Eigen::Matrix3d RotMat3(double rad)
{
    double cr = cos(rad);
    double sr = sin(rad);
    Eigen::Matrix3d Rb2a;
    Rb2a<<cr, sr, 0,  -sr, cr, 0,  0, 0, 1;
    return Rb2a;
}

//eul [R;P;Y] defined in "n": rotate "n" to obtain "b"
//result: Cb2n (from b to n) s.t., Cb2n=R3(-Y)R2(-P)R1(-R)
template<class Scalar>
static Eigen::Matrix<Scalar, 3,3 > roteu2ro(Eigen::Matrix<Scalar, 3, 1> eul)
{
    Scalar cr = cos(eul[0]); Scalar sr = sin(eul[0]);	//roll
    Scalar cp = cos(eul[1]); Scalar sp = sin(eul[1]);	//pitch
    Scalar ch = cos(eul[2]); Scalar sh = sin(eul[2]);	//heading
    Eigen::Matrix<Scalar, 3, 3> dcm;
    dcm.setZero();
    dcm(0,0) = cp * ch;
    dcm(0,1) = (sp * sr * ch) - (cr * sh);
    dcm(0,2) = (cr * sp * ch) + (sh * sr);

    dcm(1,0) = cp * sh;
    dcm(1,1) = (sr * sp * sh) + (cr * ch);
    dcm(1,2) = (cr * sp * sh) - (sr * ch);

    dcm(2,0) = -sp;
    dcm(2,1) = sr * cp;
    dcm(2,2) = cr * cp;
    return dcm;
}
#endif // SE3_FITTING_H_
