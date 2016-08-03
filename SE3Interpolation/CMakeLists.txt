SET (PROJECT_NAME SE3Interpolation)

PROJECT(${PROJECT_NAME})
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall  -O3 -march=native ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall  -O3 -march=native")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${PROJECT_SOURCE_DIR}/cmake_modules/")

find_package(Eigen3 REQUIRED)
FIND_PACKAGE(OpenCV REQUIRED)
SET(SOPHUS_INCLUDE_DIR $ENV{HOME}/ScaViSLAM/svslocal/include) # assume 1. sophus header files are put into svslocal/include/sophus directory 
# assume 2. geographiclib is installed into ScaViSLAM/svslocal direcotry

include_directories(
${EIGEN3_INCLUDE_DIR}
${SOPHUS_INCLUDE_DIR} # for both sophus and geographiclib
${OPENCV_INCLUDE_DIR}
)

LIST(APPEND LINK_LIBS
${OpenCV_LIBS}
$ENV{HOME}/ScaViSLAM/svslocal/lib/libGeographic.so
)

ADD_EXECUTABLE(${PROJECT_NAME} src/test_se3_fitting.cpp)
TARGET_LINK_LIBRARIES (${PROJECT_NAME} ${LINK_LIBS})


