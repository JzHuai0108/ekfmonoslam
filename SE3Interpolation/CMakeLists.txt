project(SE3Interpolation)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall  -O3 -march=native ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall  -O3 -march=native")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${PROJECT_SOURCE_DIR}/cmake_modules/")

find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Sophus REQUIRED)
find_package(GeographicLib REQUIRED)
add_definitions(${GeographicLib_DEFINITIONS})

include_directories(
${EIGEN3_INCLUDE_DIR}
${Sophus_INCLUDE_DIR}
${OpenCV_INCLUDE_DIRS}
${GeographicLib_INCLUDE_DIRS}
)

list(APPEND LINK_LIBS
${OpenCV_LIBS}
${GeographicLib_LIBRARIES}
)

add_executable(${PROJECT_NAME} src/test_se3_fitting.cpp)
target_link_libraries(${PROJECT_NAME} ${LINK_LIBS})


