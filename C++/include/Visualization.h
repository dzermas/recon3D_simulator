// GUARD MY HEADER FILES!
#ifndef __VISUALIZATION_H_INCLUDED__
#define __VISUALIZATION_H_INCLUDED__ 

// Local includes
#include "simulation.h"
#include "common.h"

struct color{
    short r;
    short g;
    short b;
};

class Visualization
{
private:
    // Variables
    pcl::visualization::PCLVisualizer viewer;

public:
    void EigMat2PointCloud(EigMatX&, pcl::PointCloud<pcl::PointXYZRGB>::Ptr&, color&);
    void EigVec2Point(EigVec&, pcl::PointXYZRGB&, color&);
    void VisualizePointCloud(pcl::PointCloud<pcl::PointXYZRGB>::Ptr&, std::string);
    void VisualizeCameraPose(camera&, EigMat&, EigVec&, std::string);
    void Spin();

	// Constructor / Destructor
	Visualization();
	virtual ~Visualization();
};

#endif // VISUALIZATION_H
