cmake_minimum_required(VERSION 2.8)
project(recon3D_simulator)

find_package(Eigen3)

set(EIGEN3_INCLUDE_DIR /usr/include/eigen3/Eigen)

INCLUDE_DIRECTORIES ( ${EIGEN3_INCLUDE_DIR} )


set(SOURCE_FILES ./src/main.cpp)

add_executable(simulator ${SOURCE_FILES})