//
// Created by Dimitris Zermas on 4/8/17.
//

#include "TwoImagesReconstruction.h"


/* Constructor */
TwoImagesReconstruction::TwoImagesReconstruction() {
}

/* Destructor */
TwoImagesReconstruction::~TwoImagesReconstruction(){
}
//----------------------------------------
/* SETTERS */
/* Set camera Intrinsics */
void TwoImagesReconstruction::setCameraIntrinsics(camera &cam){
    cam.fc[0] = 400; cam.fc[1] = 400;
    cam.cc[0] = 320; cam.cc[1] = 240;
    cam.kc[0] = 0; cam.kc[1] = 0; cam.kc[2] = 0; cam.kc[3] = 0; cam.kc[4] = 0;
    cam.alpha_c = 0;
    cam.dx = 640;
    cam.dy = 480;

    camera_params = cam;
}

/* Set homogeneous correspondences from the two images */
void TwoImagesReconstruction::setHomogeneousCorrespondences(EigMatX& c1, EigMatX& c2){
    c1_h = c1;
    c2_h = c2;
    data_size = c1_h.rows();
}
//----------------------------------------
/* GETTERS */
/* Get computed Essential matrix */
void TwoImagesReconstruction::getEssential(EigMat& essential){
    essential = E;
}

/* Get computed Rotation matrix */
void TwoImagesReconstruction::getRotation(EigMat& rotation){
    rotation = R;
}

/* Get computed translation vector */
void TwoImagesReconstruction::getTranslation(EigVec& translation){
    translation = t;
}

/* Get computed residual from the estimation of the essential matrix */
void TwoImagesReconstruction::getEssentialResidual(scalar& residual){
    residual = essential_residual;
}

/* Get corrected correspondences */
void TwoImagesReconstruction::getCorrectedCorrespondences(EigMatX& c1_k, EigMatX& c2_k){
    c1_k = c1_h_k;
    c2_k = c2_h_k;
}

/* Get reconstructed 3D points */
void TwoImagesReconstruction::getReconstructedPoints(EigMatX &reconstructed){
    reconstructed = r_points;
}

void TwoImagesReconstruction::getDataSize(int &size){
    size = data_size;
}
//----------------------------------------
/* FUNCTIONS */
/*  Higgins-Longuet or 8 point algorithm for estimating the Fundamental matrix
    which is the Essential Matrix here, since we have normalized the 2D
    keypoints w.r.t. the camera intrinsics.
*/
void
TwoImagesReconstruction::computeEssentialHL()
{
    EigMatX A;
    A.resize(data_size,9);
    for (size_t i = 0; i < data_size; i++){
        A(i,0) = c1_h(i,0) * c2_h(i,0);
        A(i,1) = c1_h(i,0) * c2_h(i,1);
        A(i,2) = c1_h(i,0);
        A(i,3) = c1_h(i,1) * c2_h(i,0);
        A(i,4) = c1_h(i,1) * c2_h(i,1);
        A(i,5) = c1_h(i,1);
        A(i,6) = c2_h(i,0);
        A(i,7) = c2_h(i,1);
        A(i,8) = 1;
    }

    EigVecX x(9);
    solveHomogeneousEq(A, x);

    E(0,0) = x(0);
    E(0,1) = x(1);
    E(0,2) = x(2);
    E(1,0) = x(3);
    E(1,1) = x(4);
    E(1,2) = x(5);
    E(2,0) = x(6);
    E(2,1) = x(7);
    E(2,2) = x(8);

    // SVD clean up
    Eigen::JacobiSVD<EigMatX> svd(E, Eigen::ComputeThinU | Eigen::ComputeThinV);
    EigVec D = svd.singularValues();
    EigMat U = svd.matrixU();
    EigMat V = svd.matrixV();
    D(0) = 1.0;
    D(1) = 1.0;
    D(2) = 0.0; // smallest eigenvalue is forced to be zero
    E = U * D.asDiagonal() * V.transpose();

    scalar ss = 0;
    for (size_t i = 0; i < data_size; i++){
        ss += c2_h.row(i) * E * c1_h.row(i).transpose();
    }
    essential_residual = std::abs(ss / data_size);
}

/*  Find motion parameters [R,t]
    Kanatani's methodology to compute motion parameters from the
    Fundamental/Essential matrix and the correspondences
*/
void
TwoImagesReconstruction::motionParametersKanatani()
{
    Eigen::JacobiSVD<EigMatX> svdE(E * E.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
    EigMat U = svdE.matrixU();

    t = U.col(U.cols()-1);

    //Assuming depths (Z, Z') from cameras 1 and 2 are positive, determine sign of h
    int sum = 0;
    for (size_t i = 0; i < data_size; i++){
        EigVec p1_h = c1_h.row(i).transpose();
        EigVec p2_h = c2_h.row(i).transpose();
        EigVec tmp = E * p2_h;
        sum += sgn(triple<float, EigVec>(t, p1_h, tmp));
    }
    // If the sign is negative for the t we picked, then t = -t;
    if (sum < 0) {t = -t;}

    // And compute the correct rotation R
    EigVec neg_t = -1 * t;
    EigMat K = skewSymmetric<EigMat, EigVec>(neg_t) * E;
    Eigen::JacobiSVD<EigMatX> svdK(K, Eigen::ComputeThinU | Eigen::ComputeThinV);
    EigMat Uk = svdK.matrixU();
    EigMat Vk = svdK.matrixV();

    EigVec temp = {1, 1, (Uk*Vk.transpose()).determinant()};
    R = Uk * temp.asDiagonal() * Vk.transpose();

    // Now actually see if depths are positive
    scalar Zsum = 0;
    for (size_t i = 0; i < data_size; i++){
        EigVec p1_h = c1_h.row(i).transpose();
        EigVec p2_h = c2_h.row(i).transpose();
        EigVec R_p2_h = R * p2_h;
        EigVec cross1 = cross<EigVec>(t, R_p2_h);
        EigVec cross2 = cross<EigVec>(p1_h, R_p2_h);
        Zsum += dot<scalar, EigVec>(cross1, cross2);
    }

    if (Zsum < 0) {t = -t;}
};
//----------------------------------------

/* Correct points to accurately validate the epipolar constraint, based on the first order approximation by Kanatani */
void
TwoImagesReconstruction::correctCorrespondencesKanatani()
{
	EigVec k = {0, 0, 1};
    EigMat Vv = EigMat::Identity(3,3) - k * k.transpose();
    EigMat Vu = EigMat::Identity(3,3) - k * k.transpose();

    // Statistical error of the corrected correspondences
    EigVecX J(data_size);
    // Corrected correspondences
    EigMatX c1_h_c(data_size,3);
    EigMatX c2_h_c(data_size,3);
    Dc1.resize(data_size,3);
    Dc2.resize(data_size,3);
    for (size_t i = 0; i < data_size; i++){
        EigVec v = c1_h.row(i).transpose();// v is 3x1 vector
        EigVec u = c2_h.row(i).transpose();// u is 3x1 vector
        // Compute the error of the points in the two images
        EigVec den1 = Vv * E.transpose() * v;
        scalar den1_norm = norm<scalar, EigVec>(den1);
        EigVec den2 = Vu * E * u;
        scalar den2_norm = norm<scalar, EigVec>(den2);

        scalar nom_scalar = u.transpose() * E.transpose() * v;
        scalar denominator = (den1_norm * den1_norm) + (den2_norm * den2_norm);

        EigVec nom_v_vector = Vv * E * u;
        EigVec nom_v = nom_scalar * nom_v_vector;
        Dc1.row(i).transpose() = nom_v / denominator;

        EigVec nom_u_vector = Vu * E.transpose() * v;
        EigVec nom_u = nom_scalar * nom_u_vector;
        Dc2.row(i).transpose() = nom_u / denominator;

        // Compute corrected points (data_size x 3)
        c1_h_c.row(i).transpose() = v - Dc1.row(i).transpose();
        c2_h_c.row(i).transpose() = u - Dc2.row(i).transpose();
        // Compute the residual to see whether the points correspond or not
        scalar J_nom = v.transpose() * E * u;
        EigVec J_den1 = Vv * E.transpose() * c1_h_c.row(i).transpose();
        scalar J_den1_norm = norm<scalar, EigVec>(J_den1);
        EigVec J_den2 = Vu * E * c2_h_c.row(i).transpose();
        scalar J_den2_norm = norm<scalar, EigVec>(J_den2);
        scalar J_denominator = J_den1_norm * J_den1_norm + J_den2_norm * J_den2_norm;
        J(i) = (J_nom * J_nom) / J_denominator;

    }

    // Is kept with confidence 95% if J < CritVal (Degrees of Freedom = 1)
    // For more details see Kanatani's statistical optimization book
    std::vector<int> keepers;
    for (int i = 0; i < data_size; i++){
        if (J(i) < 0.04)
            keepers.push_back(i);
    }
    // Update data_size with the number of keepers
    data_size = keepers.size();
    // Corrected correspondences that survived (keepers) after the statistical correction
    c1_h_k.resize(data_size,3);
    c2_h_k.resize(data_size,3);
    for (size_t i = 0; i < data_size; i++){
        c1_h_k.row(i) = c1_h_c.row(keepers.at(i));
        c2_h_k.row(i) = c2_h_c.row(keepers.at(i));
    }
}
//----------------------------------------

/* Triangulation based on Kanatani's method */
void
TwoImagesReconstruction::triangulationKanatani()
{
    r_points.resize(data_size,3);
    for (size_t i = 0; i < data_size; i++){
        // Compute depths of corresponding points for both cameras
        EigVec v = c1_h_k.row(i).transpose();// v is 3x1 vector
        EigVec u = c2_h_k.row(i).transpose();// u is 3x1 vector

        EigVec Ru = R * u;
        EigVec Rv = R * v;

        EigVec denom_cross = cross<EigVec>(v,Ru);
        scalar denom_cross_norm = norm<scalar,EigVec>(denom_cross);
        scalar denominator = denom_cross_norm * denom_cross_norm;

        EigVec nom_cross11 = cross<EigVec>(t,Ru);
        EigVec nom_cross12 = cross<EigVec>(v,Ru);
        EigVec nom_cross21 = cross<EigVec>(t,Rv);
        EigVec nom_cross22 = cross<EigVec>(v,Ru);
        // Depth w.r.t. camera 1
        scalar z1 = dot<scalar,EigVec>(nom_cross11, nom_cross12) / denominator;
        // Depth w.r.t. camera 1
        scalar z2 = dot<scalar,EigVec>(nom_cross21, nom_cross22) / denominator;
        EigVec n_cross = cross<EigVec>(t,v);
        EigVec n = normalize<EigVec>(n_cross);
        EigVec m = cross(n,Ru);

        // Correction of depth
        EigVec Dv = Dc1.row(i).transpose();
        EigVec Du = Dc2.row(i).transpose();
        EigVec tmp = z1*Dv - z2*R*Du;
        scalar nom_Dz = dot<scalar,EigVec>(m, tmp);
        scalar Dz = - nom_Dz / dot<scalar,EigVec>(m,v);

        r_points.row(i) = (z1 + Dz) * c1_h_k.row(i);
    }
}