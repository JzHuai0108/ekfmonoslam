# Example usage
# ./build.sh /opt/ros/kinetic/share/OpenCV-3.3.1-dev /usr/include/eigen3
OPENCV_CONFIG_PATH=$1
EIGEN3_INCLUDE_FOLDER=$2
if [[ -z $OPENCV_CONFIG_PATH ]] || [[ -z $EIGEN3_INCLUDE_FOLDER ]]; then
    echo "Usage: " $0 " OPENCV_CONFIG_PATH EIGEN3_INCLUDE_FOLDER"
    echo "Ex. " $0 " /opt/ros/kinetic/share/OpenCV-3.3.1-dev /usr/include/eigen3"
    exit -1
else
    echo "OPENCV_CONFIG_PATH " $OPENCV_CONFIG_PATH
    echo "EIGEN3_INCLUDE_FOLDER " $EIGEN3_INCLUDE_FOLDER
fi

FULL_PATH=`pwd`
INSTALL_FOLDER=$FULL_PATH/Thirdparty/slam_devel
rm -rf Thirdparty
mkdir -p $INSTALL_FOLDER

rm -rf build
mkdir -p build

cd Thirdparty
echo "Configuring and building Thirdparty/Sophus ..."
git clone https://github.com/stevenlovegrove/Sophus.git

cd Sophus
git checkout b474f05
rm -rf build
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DEIGEN3_INCLUDE_DIR:PATH=$EIGEN3_INCLUDE_FOLDER -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL_FOLDER/sophus
make -j4
make install

cd ../..
# now under Thirdparty
wget https://sourceforge.net/projects/geographiclib/files/distrib/GeographicLib-1.49.tar.gz
tar xfpz GeographicLib-1.49.tar.gz
cd GeographicLib-1.49
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$INSTALL_FOLDER/geographic
make
make install

cd ../..
# now under Thirdparty
cd ../build
cmake .. -DSophus_DIR=$INSTALL_FOLDER/sophus/share/sophus/cmake -DGeographicLib_DIR=$INSTALL_FOLDER/geographic/lib/cmake/GeographicLib
make



