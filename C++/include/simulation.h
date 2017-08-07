//
// Created by Dimitris Zermas on 4/8/17.
//

#ifndef __SIMULATION_H_INCLUDED__
#define __SIMULATION_H_INCLUDED__

//class TwoImagesReconstruction;
//class Visualization;

// C++ Libraries
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cmath>

// Eigen
#include <Eigen/Dense>

// PCL includes
#include <pcl/common/common_headers.h>
#include <pcl/common/common.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

struct camera{
    short fc[2];
    short cc[2];
    short kc[5];
    short alpha_c;
    short dx;
    short dy;
};

// Local
#include "common.h"
#include "Visualization.h"
#include "TwoImagesReconstruction.h"



#endif //SIMULATION_H
