# recon3D_simulator
Educational simulator for the reconstruction of virtual objects from a pair of images.

This example code aims to help students understand how and in what extend the various variables of a SfM scheme affect the final result.

Includes implementations of some algorithms from Kenichi Kanatani's book "Statistical Optimization for Geometric Computation".

The C++ code needs the Eigen3 and PCL-1.7 libraries.

# Eigen Installation
You may find a stable version of Eigen [here](http://eigen.tuxfamily.org/index.php?title=Main_Page), or in the terminal type:

`sudo apt install libeigen3-dev`

# PCL installation
You may find a stable version of PCL [here](http://pointclouds.org/downloads/), or in the terminal type:

`sudo apt install libpcl-1.7`

# Download
* git clone [https://github.com/dzermas/recon3D_simulator.git](https://github.com/dzermas/recon3D_simulator.git)

# Matlab version
* `cd recon3D_simulator/Matlab`
* run the `shapesSimulation.m`

# C++ version
## Compile
* `cd recon3D_simulator/C++`
* `mkdir build`
* `cd build`
* `cmake ..`
* `make`

## Run
The C++ program reads a text file with N rows and 3 columns (x,y,z) and has the following format:

x0 y0 z0

x1 y1 z1

x2 y2 z2

...

xn yn zn

To run the C++ program, build it and then from inside the build folder execute:

`./simulation ../data/<input_text_file>.txt`
