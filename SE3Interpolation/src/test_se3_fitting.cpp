
#include "se3_fitting.h"

#include <iostream>
#include <vector>
#include <fstream>

#include <Eigen/Dense> // inverse

#include <opencv2/core/core.hpp> // InspectTsukubaEulerAngles()
#include <opencv2/highgui/highgui.hpp> // InspectTsukubaEulerAngles()
#include <opencv2/imgproc.hpp> // circle()

#include <exception>
#include <GeographicLib/GravityModel.hpp> // adding gravity to simulated IMU samples
#include <GeographicLib/Geocentric.hpp>

using namespace GeographicLib;
using namespace Sophus;
using namespace std;

//input: lat and long, height is not needed
//output: Ce2n
static Matrix3d llh2dcm( Vector3d &llh)
{
    double sL = sin(llh[0]);
    double cL = cos(llh[0]);
    double sl = sin(llh[1]);
    double cl = cos(llh[1]);

    Matrix3d Ce2n;
    Ce2n<< -sL * cl, -sL * sl, cL ,  -sl, cl, 0 , -cL * cl, -cL * sl, -sL;
    return Ce2n;
}
//output rotation matrix about 3 axis for rad in radians
// if we rotate b frame about its 3 axis by rad into a frame, then =x^a*R3(rad)*x^b
static Matrix3d RotMat3(double rad)
{
    double cr = cos(rad);
    double sr = sin(rad);
    Matrix3d Rb2a;
    Rb2a<<cr, sr, 0,  -sr, cr, 0,  0, 0, 1;
    return Rb2a;
}

//eul [R;P;Y] defined in "n": rotate "n" to obtain "b"
//result: Cb2n (from b to n) s.t., Cb2n=R3(-Y)R2(-P)R1(-R)
template<class Scalar>
static Eigen::Matrix<Scalar, 3,3 > roteu2ro( Eigen::Matrix<Scalar, 3, 1> eul)
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


//convert acceleration of sensor by combined force in local world frame to accel by applied force in sensor frame
// and convert rotation velocity of sensor w.r.t world frame in sensor frame, to
// rotation velocity of sensor w.r.t inertial frame in sensor frame
// gravitySamples record gravity in e frame at each epoch
// addgnomega, do we add gravity and earth rotation into the IMU measurement? not add by default
// gravity matters much in accuracy, but earth rotation does not matter much
// fixedg, do we use constant gravity for the whole trajectory or recompute it for every point,
// note it is expensive to compute gravity every time, do not do this by default
template<class Scalar>
void ConvertToSensorFrame(const vector<Eigen::Matrix<Scalar, 4, 4> >& samplePoses, vector<Eigen::Matrix<Scalar, 10, 1> >& samples,
                          const SE3Group<Scalar> Pw2e, vector<Vector3d>& gravitySamples, bool addGnOmega=false, bool fixedG=true)
{
    typedef SE3Group<Scalar> SE3Type;
    double WIE_E=7292115e-11;    //earth's rotaion rate
    Vector3d wie2e;
    wie2e<<0.0, 0.0, WIE_E;
    int dataCount=samples.size();
    gravitySamples.resize(dataCount);    

    GravityModel grav("egm2008");
    SE3Type Ps2e;//transformation from s-frame to e-frame
    Vector3d grave;//gravity in e-frame
    double gx, gy, gz, X, Y, Z;

    //compute gravity in e frame
    X=Pw2e.translation()[0];
    Y=Pw2e.translation()[1];
    Z=Pw2e.translation()[2];
    grav.W(X,Y,Z, gx, gy, gz);// gravity in e-frame
    grave<<gx,gy,gz;
    if(!addGnOmega){
          grave<<0,0,0;
          wie2e<<0,0,0;
    }

    for( int j=0; j<dataCount; ++j)
    {
        Ps2e=Pw2e*SE3Type(samplePoses[j]);
        if(addGnOmega&&(!fixedG)){
            //compute gravity in e frame
            X=Ps2e.translation()[0];
            Y=Ps2e.translation()[1];
            Z=Ps2e.translation()[2];
            grav.W(X,Y,Z, gx, gy, gz);// gravity in e-frame
            grave<<gx,gy,gz;
            //the following code snippet tests gravity computation in another way gets the same result
            /*cout<<"W:"<< gx<< ","<< gy<<", "<<gz<<endl;
            double lat, lon, h;
            earth.Reverse(X,Y,Z,lat, lon, h);

            grav.Gravity(lat, lon, h, gx, gy, gz);//gravity in n-frame
            cout << "gravity:"<< gx << " " << gy << " " << gz << "\n";

            Vector3d inillh, gravn;
            inillh<<lat*M_PI/180, lon*M_PI/180, h;
            cout<<inillh<<endl;
            gravn<< gy, gx, -gz;
            SO3Group<double> Re2n0=llh2dcm(inillh);
            Vector3d grave0=Re2n0.inverse()*gravn;
            cout<<"grave:"<< grave0<<endl;*/
        }
        gravitySamples[j]=grave;

#if 0 // use average acceleration. Empricially, this does not make much difference as true epoch acceleration
        if(j==0)
            samples[j].segment(1,3)=samplePoses[j].block(0,0,3,3).transpose()*samples[j].segment(1,3)-
                Ps2e.rotationMatrix().inverse()*grave;
        else
        {
            Scalar interval=(samples[j][0]-samples[j-1][0]);
            samples[j].segment(1,3)=samplePoses[j].block(0,0,3,3).transpose()*
                ((samples[j].tail(3)-samples[j-1].tail(3))/interval)-
                Ps2e.rotationMatrix().inverse()*grave;
        }
#else
        samples[j].segment(1,3)=samplePoses[j].block(0,0,3,3).transpose()*samples[j].segment(1,3)-
            Ps2e.rotationMatrix().inverse()*grave;
#endif
        samples[j].segment(4,3)=samples[j].segment(4,3)+Ps2e.rotationMatrix().inverse()*wie2e;
    }
}


// this function is written by Jianzhu Huai, Dec 19 2014 to construct
// cubic B splines given a series of control poses as done in Lovegrove 2013
// then angular rate, $\omega_{ws}^s$ and acceleration $a_s^w$ are interpolated
// assume the input data poses are almost uniformly distributed in time

template<class Scalar>
void SimulateIMU(int testCase=0)
{
    typedef SO3Group<Scalar> SO3Type;
    typedef SE3Group<Scalar> SE3Type;
    typedef typename SE3Group<Scalar>::Point Point;
    typedef Eigen::Matrix<Scalar, 10, 1> TAO; // time, acceleration, omega/angular rate, and velocity

    cerr<< "my test on se3 interpolation with B splines"<<endl;
    cerr<< "setting sampling parameters"<<endl;
    string poseFile, timeFile;
    string outputFile, debugPoseFile, initFile;
    Scalar outputFreq;
    bool fixedG, addGnOmega;
    Scalar lat(0.0), lon(0.0), h(0.0);
    SO3Type Rw2n0; // w-frame is the first IMU frame, n0 frame is the n-frame at that epoch
    Eigen::Matrix<Scalar, 3,3 > Rs2c;//the rotation matrix from IMU sensor frame to camera frame
    Point tsinc;//the coordinate of sensor frame origin in camera frame
    SE3Type Ts2c, Tc2s;
    vector<SE3Type> q02n; //q_0^w, q_1^w, ..., q_n^w; N=n+1 poses, sensor frame to world frame transformations

    vector<Scalar > times; //timestamps for N poses
    string tempStr;
    if(testCase==0){ //tailored for KITTI odometry dataset color seq 00
        poseFile="/media/jhuai/Mag/KITTISeq00/00.txt";//stores the transform Tc2w, camera to world
        timeFile="/media/jhuai/Mag/KITTISeq00/times.txt";
        outputFreq=100; // the higher frequency, the better reconstructed trajectory,
        // gravity models also makes a difference for KITTI color seq 00
        outputFile="/media/jhuai/Mag/KITTISeq00/SimulateIMUTest/samples.txt";
        debugPoseFile="/media/jhuai/Mag/KITTISeq00/SimulateIMUTest/poses.txt";
        initFile="/media/jhuai/Mag/KITTISeq00/SimulateIMUTest/initCond.txt";
        fixedG=false;
        addGnOmega=true;

        lat = 8.389733;
        lon = 48.984678;
        h = 110.0; // Karlsruhe University, roughly estimated from
        // Geiger, Andreas, et al. "Vision meets robotics: The KITTI dataset."
        // The International Journal of Robotics Research (2013): 0278364913491297. and OpenStreetMap
        Rw2n0=RotMat3(-15*M_PI/180);

        Rs2c<<0,1,0,0,0,1,1,0,0;
        tsinc=Point(0.05, 0.1, -0.1);//the coordinate of sensor frame origin in left P0 camera frame

        cerr<< "read in data poses"<<endl;

        Ts2c=SE3Type(SO3Type(Rs2c), tsinc);
        Tc2s=Ts2c.inverse();

        ifstream dataptr(poseFile.c_str());
        assert(!dataptr.fail());

        Eigen::Matrix<Scalar,4,4> transMat;
        Scalar precursor=0;
        int lineNum=0;
        while(!dataptr.eof()){
            dataptr>>precursor;
            if(dataptr.fail())
                break;
            transMat.setZero();
            transMat.data()[0]=precursor;
            for (int j=1; j<12; ++j)
                dataptr>>transMat.data()[j];
            ++lineNum;
            q02n.push_back(Tc2s*SE3Type(transMat.transpose())*Ts2c); //Ts2s0
            //getline(dataptr, tempStr);
        }
        dataptr.close();
        cout<<"Num of lines:"<<lineNum<<endl;
        ifstream timeptr(timeFile.c_str());
        assert(!timeptr.fail());
        times.resize(lineNum);

        int counter=0;
        while(!timeptr.eof()){
            timeptr>>precursor;
            if(timeptr.fail())
                break;
            times[counter]=precursor;
            ++counter;
            //getline(dataptr, tempStr);
        }
        timeptr.close();
        assert(lineNum==counter);
    }
    else if(testCase==1)
    {
        //tailored for Tsukuba 3G stereo dataset
        poseFile="/media/jhuai/Mag/NewTsukubaStereoDataset/groundtruth/camera_track.txt";
        //stores the transform Tc2w, stereo camera central frame to world, X,Y,Z and A,B,C of Euler Angles
        Scalar inputFreq=30;
        outputFreq=100; // the higher frequency, the better reconstructed trajectory,
        outputFile="/media/jhuai/Mag/NewTsukubaStereoDataset/SimulateIMUTest/samples.txt";
        debugPoseFile="/media/jhuai/Mag/NewTsukubaStereoDataset/SimulateIMUTest/poses.txt";
        initFile="/media/jhuai/Mag/NewTsukubaStereoDataset/SimulateIMUTest/initCond.txt";
        fixedG=false;
        addGnOmega=true;

        lat = 40.003227423;
        lon = -83.042986648;
        h = 212.6010; //Ohio State University
        Eigen::Matrix<Scalar, 3,3> tempRot;
        tempRot<<1,0,0,0,1,0,0,0,1;
        Rw2n0=tempRot;
        Rs2c<<0,1,0,0,0,-1,-1,0,0;
        tsinc=Point(0.0, -0.06, 0.02);//the coordinate of sensor frame origin in central camera frame

        cerr<< "read in data poses"<<endl;
        Ts2c=SE3Type(SO3Type(Rs2c), tsinc);
        Tc2s=Ts2c.inverse();

        ifstream dataptr(poseFile.c_str());
        assert(!dataptr.fail());

        Eigen::Matrix<Scalar,4,4> transMat;
        Point eulc2w;
        int lineNum=0;
        while(getline(dataptr,tempStr))
        {
            std::stringstream   linestream(tempStr);
            std::string        value;
            int valueNum=0;
            while(getline(linestream,value,','))
            {
                if(valueNum<3)
                    transMat.data()[valueNum+12]=atof(value.c_str())/100;//convert from cm to m
                else
                    eulc2w[valueNum-3]=atof(value.c_str())*M_PI/180;
                ++valueNum;
            }
            transMat.topLeftCorner(3,3)=roteu2ro(eulc2w);
            ++lineNum;
            q02n.push_back(Tc2s*SE3Type(transMat)*Ts2c); // this is Ts2s0
        }
        dataptr.close();
        cout<<"Num of lines:"<<lineNum<<endl;

        times.resize(lineNum);
        for(int j=0; j<lineNum;++j)
        {
            times[j]=j/inputFreq;
        }
    }
    vector<Eigen::Matrix<Scalar, 4,4 > > samplePoses;
    vector<TAO> samples;
    InterpolateIMUData(q02n, times, outputFreq, samplePoses,  samples);

    //add gravity and earth rotation effect, optional
    Geocentric earth(Constants::WGS84_a(), Constants::WGS84_f());
    // Alternatively: const Geocentric& earth = Geocentric::WGS84();

    double X, Y, Z;
    earth.Forward(lat, lon, h, X, Y, Z);
    Vector3d inillh;
    inillh<<lat*M_PI/180, lon*M_PI/180, h;
    SO3Type Re2n0=llh2dcm(inillh);

    Point twine(X,Y,Z);
    SE3Type Pw2e(Re2n0.inverse()*Rw2n0, twine);

    ofstream poseptr(debugPoseFile);

    ofstream sampleptr(outputFile);
    sampleptr<<"%% timestamp(sec), $a_{is}^s$(m/s^2), $\\Omega_{is}^s$(rad/sec), $v_{ws}^w$(m/s), $gravity^e$."<<endl;
    ofstream initptr(initFile);
    initptr<<"initial lat(deg), lon(deg), height(m) in WGS84 CRS"<<endl;
    initptr.precision(10);
    initptr<<Point(lat,lon,h)<<endl;
    initptr<<"initial rotation matrix from IMU sensor frame to NED frame"<<endl;
    initptr<<Rw2n0.matrix()<<endl;
    initptr<<"Ts2c transformation matrix from IMU sensor frame to camera frame used in POSE INPUT file(row major)"<<endl;
    initptr<<Ts2c.matrix().block(0,0,3,4)<<endl;
    initptr.close();
    vector<Vector3d> gravitySamples;
    ConvertToSensorFrame(samplePoses,samples, Pw2e, gravitySamples, addGnOmega, fixedG);
    int dataCount=samplePoses.size();
    for (int i=0; i<dataCount;++i){
        poseptr<<samples[i][0]<<" "<<samplePoses[i].row(0)<<" "<<samplePoses[i].row(1)<<" ";
        poseptr<<samplePoses[i].row(2)<<endl;
        sampleptr<<samples[i].transpose()<<" "<<gravitySamples[i].transpose()<<endl;
    }
    sampleptr.close();
    Eigen::Matrix<Scalar, 4,4> tempTrans=q02n.rbegin()->matrix();
    poseptr<<*(times.rbegin())<<" "<<tempTrans.row(0)<<" "<<tempTrans.row(1)<<" "<<tempTrans.row(2)<<endl;
    poseptr.close();
}
// Example of using the GeographicLib::GravityModel class
int TestGravity() {
    try {

        GravityModel grav("egm2008");

        double lat = 27.99, lon = 86.93, h = 8820; // Mt Everest
        double X, Y, Z;
        Geocentric earth(Constants::WGS84_a(), Constants::WGS84_f());
        earth.Forward(lat, lon, h, X, Y, Z);

        double gx, gy, gz;
        grav.W(X,Y,Z, gx, gy, gz);//gravity model, W=V+Phi
        cout << "W:"<< gx << " " << gy << " " << gz << "\n";
        grav.V(X,Y,Z, gx, gy, gz);//gravitational model egm2008 or else
        Vector3d grave;
        grave<<gx, gy, gz;
        cout <<"V:"<< gx << " " << gy << " " << gz << "\n";
        grav.Gravity(lat,lon, h, gx, gy, gz);// the same as W
        cout <<"gravity:"<< gx << " " << gy << " " << gz << "\n";
        Vector3d inillh, gravn;
        inillh<<lat*M_PI/180, lon*M_PI/180, h;
        gravn<< gy, gx, -gz;
        SO3Group<double> Re2n0=llh2dcm(inillh);
        cout<< "gravn:"<<gravn<<endl;
        Vector3d grave0=Re2n0.inverse()*gravn;
        cout<<"grave:"<< grave0<<endl;
    }
    catch (const exception& e) {
        cerr << "Caught exception: " << e.what() << "\n";
        return 1;
    }
    try {
        Geocentric earth(Constants::WGS84_a(), Constants::WGS84_f());
        // Alternatively: const Geocentric& earth = Geocentric::WGS84();
        {
            // Sample forward calculation
            double lat = 27.99, lon = 86.93, h = 8820; // Mt Everest
            double X, Y, Z;
            earth.Forward(lat, lon, h, X, Y, Z);
            cout << floor(X / 1000 + 0.5) << " "
                 << floor(Y / 1000 + 0.5) << " "
                 << floor(Z / 1000 + 0.5) << "\n";
        }
        {
            // Sample reverse calculation
            double X = 302e3, Y = 5636e3, Z = 2980e3;
            double lat, lon, h;
            earth.Reverse(X, Y, Z, lat, lon, h);
            cout << lat << " " << lon << " " << h << "\n";
        }
    }
    catch (const exception& e) {
        cerr << "Caught exception: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
void InspectTsukubaEulerAngles()
// the Tsukuba dataset gives unclear euler angles, this function displays trajectory
// of some point over frames, verifies the following:
//(1) Tsukuba camera coordinate system is defined as right, up and back for x,y,z
// when the we are viewing from behind of the camera
//(2) In camera_track.txt, T_{c_i}^{c_0} is given in each line, such that t_{c_i}^{c_0}=[X,Y,Z]',
// and R_{c_i}^{c_0}=R3(-C)R2(-B)R1(-A), and T_{c_i}^{c_0}=[R_{c_i}^{c_0}, t_{c_i}^{c_0}]
//(3) Its disparity map for occluded points are invalid
{
    using namespace cv;
    typedef SE3Group<double> SE3Type;
    //Tsukuba 3G stereo dataset
    string poseFile="/media/jhuai/Mag/NewTsukubaStereoDataset/groundtruth/camera_track.txt";
    //stores the transform Tc2w, stereo camera central frame to world, X,Y,Z and A,B,C of Euler Angles
    string leftImagePath="/media/jhuai/Mag/NewTsukubaStereoDataset/illumination/daylight/left/tsukuba_daylight_L_";
    char leftImageName[300]={'\0'};
    string rightImagePath="/media/jhuai/Mag/NewTsukubaStereoDataset/illumination/daylight/right/tsukuba_daylight_R_";
    char rightImageName[300]={'\0'};

    string leftImageDepth="/media/jhuai/Mag/NewTsukubaStereoDataset/groundtruth/depth_maps/left/tsukuba_depth_L_00001.xml";
    string rightImageDepth="/media/jhuai/Mag/NewTsukubaStereoDataset/groundtruth/depth_maps/right/tsukuba_depth_R_00001.xml";
    string leftImageDisp="/media/jhuai/Mag/NewTsukubaStereoDataset/groundtruth/disparity_maps/left/tsukuba_disparity_L_00001.png";
    string rightImageDisp="/media/jhuai/Mag/NewTsukubaStereoDataset/groundtruth/disparity_maps/right/tsukuba_disparity_R_00001.png";
    Vector3d point, imageObs;
    Matrix3d K;
    K<<615, 0, 320, 0, 615, 240, 0, 0, 1;
    imageObs<<600, 100, 1;
    point=K.inverse()*imageObs;
    cv::Mat leftDepth, rightDepth;
    cv::FileStorage fs;
    fs.open(leftImageDepth, cv::FileStorage::READ);
    if (!fs.isOpened())
    {
        cerr << "Failed to open " << leftImageDepth << endl;
        return;
    }
    fs["depth"] >> leftDepth;
    fs.release();
    fs.open(rightImageDepth, cv::FileStorage::READ);
    if (!fs.isOpened())
    {
        cerr << "Failed to open " << rightImageDepth << endl;
        return;
    }
    fs["depth"] >> rightDepth;
    fs.release();

    cout<<"left depth:"<<(leftDepth.at<float>(imageObs[1], imageObs[0])/100)<<endl;
    point*=(leftDepth.at<float>(imageObs[1], imageObs[0])/100);
    point[0]-=0.05;
    point[1]=-point[1];
    point[2]=-point[2];

    Mat leftDisp=cv::imread(leftImageDisp);
    Mat rightDisp=cv::imread(rightImageDisp);
    int dispL=leftDisp.at<uchar>(imageObs[1], imageObs[0]);
    cout<<"dispL:"<<dispL<<endl;
    cout<<"right depth:"<<rightDepth.at<float>(imageObs[1], imageObs[0]-dispL)/100<<endl;
    int dispR=rightDisp.at<uchar>(imageObs[1], imageObs[0]-dispL);
    cout<<"dispR:"<<dispR<<endl;
    sprintf(leftImageName, "%s%05d.png", leftImagePath.c_str(), 1);
    cv::Mat image=cv::imread(leftImageName);
//occlusion points may not correspond well according to their disparity
    cv::circle(image, Point(imageObs[0], imageObs[1]), 3, cv::Scalar(0,0,255), 2, 8, 0);
    cv::imshow("Point left", image);
    sprintf(rightImageName, "%s%05d.png", rightImagePath.c_str(), 1);
    image=cv::imread(rightImageName);

    cv::circle(image, cv::Point(imageObs[0]-dispL, imageObs[1]),
            3, cv::Scalar(0,0,255), 2, 8, 0);
    cv::imshow("Point right", image);
    cv::waitKey(0);
    cv::destroyAllWindows();
    std::vector<cv::Point2f> imagePointTraj;

    vector<SE3Type> q02n; //T_{c_i}^{c_0}, i=0,..., n, N=n+1

    Eigen::Matrix<double,4,4> transMat;
    transMat<<1,0,0,0.05,
              0,-1,0,0,
              0,0,-1,0,
              0,0,0,1;
    SE3Type Tc2l(transMat); //transformation from center frame to left frame
    cerr<< "read in data poses"<<endl;
    ifstream dataptr(poseFile.c_str());
    assert(!dataptr.fail());
    string tempStr;

    Eigen::Matrix<double, 3,1 > eulc2w;
    int lineNum=0;
    while(getline(dataptr,tempStr))
    {
        std::stringstream   linestream(tempStr);
        std::string        value;

        int valueNum=0;
        while(getline(linestream,value,','))
        {
            if(valueNum<3)
                transMat.data()[valueNum+12]=atof(value.c_str())/100;//convert from cm to m
            else
                eulc2w[valueNum-3]=atof(value.c_str())*M_PI/180;
            ++valueNum;
        }
        transMat.topLeftCorner(3,3)=roteu2ro(eulc2w);

        ++lineNum;
        q02n.push_back(SE3Type(transMat)); // this is Ts2s0
    }
    dataptr.close();
    cout<<"Num of lines:"<<lineNum<<endl;
    imagePointTraj.resize(lineNum);
    Vector3d tempPt;
    for( int j=0; j<lineNum; ++j)
    {
        tempPt=Tc2l*q02n[j].inverse()*point;
        if(tempPt[2]<0 || abs(tempPt[2])<1e-6)
        {
            imagePointTraj[j].x=-10;
            imagePointTraj[j].y=-10;
        }else{
        tempPt/=tempPt[2];
        tempPt=K*tempPt;
        imagePointTraj[j].x=tempPt[0];
        imagePointTraj[j].y=tempPt[1];
        }
    }
    for( int j=1; j<lineNum+1; ++j)
    {
        sprintf(leftImageName, "%s%05d.png", leftImagePath.c_str(), j);
        image=cv::imread(leftImageName);

        cv::circle(image, imagePointTraj[j-1], 3, cv::Scalar(0,0,255), 2, 8, 0);
        cv::imshow("Point Trajectory", image);
        cv::waitKey(33);
    }
}


int main() {
    //TestGravity();
    //InspectTsukubaEulerAngles();
    SimulateIMU<double>(0);
    return 0;
}
