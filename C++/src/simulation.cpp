#include <simulation.h>

using std::cout;
using std::endl;
using std::vector;

using namespace Eigen;

int main(int argc, char *argv[]){

    if (argc != 2){
        cout << "Try: ./simulation ../data/<data_file>.txt" << endl;
        exit(1);
    }
    // Initialize the random number generator
    srand(time(NULL));

    const scalar CAM_T_X = 2.0;
    const scalar CAM_T_Y = 1.0;
    const scalar CAM_T_Z = 0.2;
    const scalar NOISE = 0.00;

	// Camera Intrinsics
	camera cam; // Initialized in the simulation.h file

	/*Create random motion parameters*/
	// Rotations
	// reference R1
    EigMat C1_R_C1 = EigMat::Identity(3,3);
	// reference R2
    EigMat C2_R_C1 = generateRandomRotation<EigMat, scalar>();

	// Translations
	// reference t1
    EigVec C1_P_C1(0, 0, 0);
	// reference t2
    EigVec C2_P_C1;
	C2_P_C1 << CAM_T_X*rand01<scalar>(), CAM_T_Y*rand01<scalar>(), CAM_T_Z*rand01<scalar>();

    // Reference Essential matrix (can be used as reference for the user)
    EigMat E_ref;
    E_ref = skewSymmetric<EigMat, EigVec>(C2_P_C1) * C2_R_C1;
    //printMatrix(E_ref);

    // Read the 3D points
    EigMatX C1_P; //Reference 3D points (euclidean values)
    readData(argv[1], C1_P);

    // Initialize all related data matrices
    MatrixXf C1_pix(C1_P.rows(),2); // Pixel (pixel values)
    EigMatX C1_h(C1_P.rows(),3); // Homogeneous (normalized euclidean values)

    EigMatX C2_P(C1_P.rows(),3); //3D points w.r.t. the second camera (euclidean values)
    MatrixXf C2_pix(C1_P.rows(),2); // Pixel (pixel values)
    EigMatX C2_h(C1_P.rows(),3); // Homogeneous (normalized euclidean values)

    // Generate points on the image planes
    // Try to see how much is considered "too much noise" :)
    for (short i = 0; i < C1_P.rows(); i++){
        EigVec euclidean = C1_P.row(i);

        // Transform 3D point to homogeneous coordinates
        EigVec homogeneous;
        homogeneous(0) = euclidean(0) / euclidean(2);
        homogeneous(1) = euclidean(1) / euclidean(2);
        homogeneous(2) = euclidean(2) / euclidean(2);
        // Transform homogeneous coordinates to pixel coordinates of image 1
        Vector2f pixels;
        pixels(0) = homogeneous(0) * cam.fc[0] + cam.cc[0];
        pixels(1) = homogeneous(1) * cam.fc[1] + cam.cc[1];
        // Gaussian noise in measurements of image 1
        EigVec g_noise_1 = {NOISE * rand01<scalar>(), NOISE * rand01<scalar>(), 0.0};
        // Add noise to the homogeneous coordinates
        EigVec homogeneous_noisy = homogeneous + g_noise_1;
        // Store the noisy homogeneous coordinates
        C1_h.row(i) = homogeneous_noisy;
        // Store the noisy pixel values of image 1
        C1_pix.row(i) = pixels;
        // Transform values of the coordinate system of image 1 to the coordinate system of image 2
        // Actual 3D points w.r.t. image 2
        EigVec C2_euclidean = C2_R_C1.transpose() * euclidean - C2_R_C1 * C2_P_C1;
        // Transformed into homogeneous coordinates
        EigVec C2_homogeneous;
        C2_homogeneous(0) = C2_euclidean(0) / C2_euclidean(2);
        C2_homogeneous(1) = C2_euclidean(1) / C2_euclidean(2);
        C2_homogeneous(2) = C2_euclidean(2) / C2_euclidean(2);
        // Transformed into pixels of the image 2
        Vector2f C2_pixels;
        C2_pixels(0) = C2_homogeneous(0) * cam.fc[0] + cam.cc[0];
        C2_pixels(1) = C2_homogeneous(1) * cam.fc[1] + cam.cc[1];

        C2_P.row(i) = C2_euclidean;

        C2_h.row(i) = C2_homogeneous;

        C2_pix.row(i) = C2_pixels;
    }

    TwoImagesReconstruction tir;
    tir.setCameraIntrinsics(cam);
    tir.setHomogeneousCorrespondences(C1_h, C2_h);

    tir.computeEssentialHL();
    tir.motionParametersKanatani();
    tir.correctCorrespondencesKanatani();
    tir.triangulationKanatani();
    EigMatX reconstructed_points;
    tir.getReconstructedPoints(reconstructed_points);
    //printMatrix(reconstructed_points);

    /* Plot 3D */
    Visualization viewer;
    // Visualize camera 1
    viewer.VisualizeCameraPose(cam, C1_R_C1, C1_P_C1, "C1");
    // Visualize camera 2 reference
    viewer.VisualizeCameraPose(cam, C2_R_C1, C2_P_C1, "C2_Ref");
    // Visualize camera 2
    EigMat R; tir.getRotation(R);
    EigVec t; tir.getTranslation(t);
    viewer.VisualizeCameraPose(cam, R, t, "C2");

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr ref_cloud (new pcl::PointCloud<pcl::PointXYZRGB>);
    color rgb_ref = {.r = 255, .g = 0, .b = 0}; // initialize color struct to red
    viewer.EigMat2PointCloud(C1_P, ref_cloud, rgb_ref);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr rec_cloud (new pcl::PointCloud<pcl::PointXYZRGB>);
    color rgb_rec = {.r = 0, .g = 0, .b = 255}; // initialize color struct to blue
    viewer.EigMat2PointCloud(reconstructed_points, rec_cloud, rgb_rec);

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr final_cloud (new pcl::PointCloud<pcl::PointXYZRGB>);
    *final_cloud += *ref_cloud;
    *final_cloud += *rec_cloud;
    viewer.VisualizePointCloud(final_cloud, "3D Viewer");

    viewer.Spin();

	return 1;
}
