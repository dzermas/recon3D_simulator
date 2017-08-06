# recon3D_simulator
Educational simulator for the reconstruction of virtual objects from a pair of images.

This example code aims to help students understand how and in what extend the various variables of a SfM scheme affect the final result.

Includes implementations of some algorithms from Kenichi Kanatani's book "Statistical Optimization for Geometric Computation".

The C++ code needs the Eigen and VTK libraries.

# Eigen Installation
You may find a stable version of Eigen [here](http://eigen.tuxfamily.org/index.php?title=Main_Page).

# VTK installation
You may find a stable version of VTK [here](https://github.com/Kitware/VTK/tree/release-6.3) and instructions in building it [here](http://www.vtk.org/Wiki/VTK/Building/Linux).

# Download
* git clone [https://github.com/dzermas/recon3D_simulator.git](https://github.com/dzermas/recon3D_simulator.git)

# Matlab version
* cd recon3D_simulator/Matlab
* run the "shapesSimulation.m"

# C++ version
## Compile
* cd recon3D_simulator/C++
* mkdir build
* cd build
* cmake ..
* make

## Run
./simulation ../data/<input_text_file>.txt
