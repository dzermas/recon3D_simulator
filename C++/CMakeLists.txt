cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

## Project name
project(Simulation3D)

## C++ Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

## Static Variables
set(Local_INCLUDE_DIRS include)
set(SOURCE src/simulation.cpp src/TwoImagesReconstruction.cpp src/Visualization.cpp)
#set(VTK_DIR )

## Packages
find_package(Eigen3 REQUIRED)
find_package(PCL 1.7 REQUIRED)
find_package(VTK REQUIRED)

## External and Local includes
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${PCL_INCLUDE_DIRS})
include_directories(${Local_INCLUDE_DIRS})

## PCL
add_definitions(${PCL_DEFINITIONS})
link_directories(${PCL_LIBRARY_DIRS})

## Executables
add_executable (simulation ${SOURCE})

## External Libraries
target_link_libraries (simulation ${PCL_LIBRARIES} ${VTK_LIBRARIES})
