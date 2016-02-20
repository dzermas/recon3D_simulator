//
// Created by dzermas on 2/17/16.
//
#include <iostream>
#include <Eigen>

using namespace Eigen;
using namespace std;

// Class
struct point{
    Vector3d euclidean;
    Vector3d homogeneous;
    Vector3d bearing;
    Vector2d pixel;
    double   depth;
};

struct cameraParameters{
    Vector2i fc;
    Vector2i cc;
    int      dx;
    int      dy;
};

// Function declaration
point            initializePoint(int numPoints, double depthRange, cameraParameters camera);
cameraParameters initializeCamera();


int
main(int argc, char** argv) {
    // Initialize parameters
    cameraParameters camera;
    int numPoints = 10;
    double depthRange = 1.0;

    VectorXd   c_d(numPoints);      // point depth
    Matrix3Xcd c1_h(numPoints, 3);  // homogeneous coordinates for camera1
    Matrix3Xcd c1_b(numPoints, 3);  // bearings for camera1
    Matrix3Xcd c1_P(numPoints, 3);  // 3D point w.r.t. camera1
    Matrix3Xcd c1_pix(numPoints, 3);// pixel coordinates in camera1

    Matrix3Xcd c2_h(numPoints, 3);  // homogeneous coordinates for camera2
    Matrix3Xcd c2_b(numPoints, 3);  // bearings for camera2
    Matrix3Xcd c2_P(numPoints, 3);  // 3D point w.r.t. camera2
    Matrix3Xcd c2_pix(numPoints, 3);// pixel coordinates in camera2

    // Initialize camera parameters
    camera = initializeCamera();    // TODO: Read a config file

    // Populate the corresponding points between the two cameras
    for (int it = 0; it < numPoints; it++) {

        point c1; point c2;
        initializePoint(numPoints, depthRange, camera, c1, c2);

        c_d.row(it) << c1.depth;
        c1_h.row(it) << c1.homogeneous;
        c1_b.row(it) << c1.bearing;
        c1_P.row(it) << c1.euclidean;
        c1_pix.row(it) << c1.pixel;

    }
    std::cout << camera.fc(0) << std::endl;
    return 0;
}

cameraParameters initializeCamera() {
    cameraParameters camera;
    camera.fc(0) = 400; camera.fc(1) = 400;
    camera.cc(0) = 320; camera.cc(1) = 240;
    camera.dx = 640;
    camera.dy = 480;
    return camera;
}

void
initializePoint(int numPoints, double depthRange, cameraParameters camera, point c1, point c2){
    double randPixelX = round(rand() * camera.dx);
    double randPixelY = round(rand() * camera.dy);
    double randDepth  = 5 + round(rand() * depthRange);

    c1.pixel(0) = randPixelX; c1.pixel(1) = randPixelY;
    c1.homogeneous = (c1.pixel - camera.cc) / camera.fc  ;
    c1.euclidean


}