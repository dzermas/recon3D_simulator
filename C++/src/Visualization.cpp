//
// Created by Dimitris Zermas on 6/8/17.
//

#include "Visualization.h"

/* Constructor */
Visualization::Visualization(){
    viewer.setBackgroundColor (0, 0, 0);
    viewer.addCoordinateSystem (1.0);
    viewer.initCameraParameters ();
}

/* Destructor */
Visualization::~Visualization(){
}

void Visualization::EigMat2PointCloud(EigMatX& pts, pcl::PointCloud<pcl::PointXYZRGB>::Ptr& point_cloud, color& rgb){

    int data_size = pts.rows();
    // Create point cloud and fill it with our data
    point_cloud->resize(data_size);
    for (int i = 0; i < data_size; i++){
        point_cloud->at(i).x = pts(i,0);
        point_cloud->at(i).y = pts(i,1);
        point_cloud->at(i).z = pts(i,2);
        point_cloud->at(i).r = rgb.r;
        point_cloud->at(i).g = rgb.g;
        point_cloud->at(i).b = rgb.b;
    }
}

void Visualization::EigVec2Point(EigVec& vec, pcl::PointXYZRGB& pt, color& rgb){
    pt.x = vec(0);
    pt.y = vec(1);
    pt.z = vec(2);
    pt.r = rgb.r;
    pt.g = rgb.g;
    pt.b = rgb.b;
}

/* Plot 3D points */
void Visualization::VisualizePointCloud(pcl::PointCloud<pcl::PointXYZRGB>::Ptr& point_cloud, std::string cloud_name){
    pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(point_cloud);
    viewer.addPointCloud<pcl::PointXYZRGB> (point_cloud, rgb, cloud_name);
    viewer.setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, cloud_name);
}

void Visualization::VisualizeCameraPose(camera& cam, EigMat& R, EigVec& t, std::string id){
    scalar fc = cam.fc[0];
    scalar w = cam.cc[0] / fc;
    scalar h = cam.cc[1] / fc;
    scalar f = 1;

    EigMatX cam_face(3,4), Rt(3,4), face_edges(4,4);
    Rt << R, t;
    face_edges << -w,  w,  w, -w,
                   h,  h, -h, -h,
                   f,  f,  f,  f,
                   1,  1,  1,  1;

    cam_face = Rt * face_edges;

    // Assign random color ro cloud
    color rgb;
    rgb.r = 0;
    rgb.g = 1;
    rgb.b = 0;
    // Turn Eigen points to PCL points
    pcl::PointXYZRGB f1, f2, f3, f4, cam_base;
    EigVec tmp = cam_face.block(0,0,3,4).col(0);
    EigVec2Point(tmp, f1, rgb);
    tmp = cam_face.block(0,0,3,4).col(1);
    EigVec2Point(tmp, f2, rgb);
    tmp = cam_face.block(0,0,3,4).col(2);
    EigVec2Point(tmp, f3, rgb);
    tmp = cam_face.block(0,0,3,4).col(3);
    EigVec2Point(tmp, f4, rgb);

    EigVec2Point(t, cam_base, rgb);
    // Create camera face
    viewer.addLine(f1, f2, rgb.r, rgb.g, rgb.b, "f1" + id);
    viewer.addLine(f2, f3, rgb.r, rgb.g, rgb.b, "f2" + id);
    viewer.addLine(f3, f4, rgb.r, rgb.g, rgb.b, "f3" + id);
    viewer.addLine(f4, f1, rgb.r, rgb.g, rgb.b, "f4" + id);
    // Create camera base
    viewer.addLine(f1, cam_base, rgb.r, rgb.g, rgb.b, "base1" + id);
    viewer.addLine(f2, cam_base, rgb.r, rgb.g, rgb.b, "base2" + id);
    viewer.addLine(f3, cam_base, rgb.r, rgb.g, rgb.b, "base3" + id);
    viewer.addLine(f4, cam_base, rgb.r, rgb.g, rgb.b, "base4" + id);
    // Add camera's name
    cam_base.z -= 0.05; // move name tag a bit for better visualization. Also, when on {0,0,0} the addText3D does not work
    viewer.addText3D(id, cam_base, 0.2, 1, 1, 1, id);
}

void Visualization::Spin(){
    while (!viewer.wasStopped ())
    {
        viewer.spinOnce (100);
    }

    std::cout << "Viewer terminated" << std::endl;
    exit(1);

}