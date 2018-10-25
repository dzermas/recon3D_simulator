// GUARD MY HEADER FILES!
#ifndef __TWOIMAGESRECONSTRUCTION_H_INCLUDED__
#define __TWOIMAGESRECONSTRUCTION_H_INCLUDED__

#include "simulation.h"
#include "common.h"

/* Class implementing the 3D reconstruction from 2 images algorithm */
class TwoImagesReconstruction
{
private:
    // Variables
    camera camera_params;

    EigMatX c1_h; // homogeneous coordinates of image 1
    EigMatX c2_h; // homogeneous coordinates of image 2

    EigMatX c1_h_k; // corrected homogeneous coordinates of image 1
    EigMatX c2_h_k; // corrected homogeneous coordinates of image 2

    EigMatX Dc1; // Correction of the coordinates in image 1
    EigMatX Dc2; // Correction of the coordinates in image 2

    EigMat E; // Computed Essential matrix
    EigMat R; // Computed rotation matrix
    EigVec t; // Computed translation vector

    EigMatX r_points; // Reconstructed 3D points

    scalar essential_residual; // Residual from the estimation of the essential matrix

    int data_size;
public:
	// Set functions
	void setCameraIntrinsics(camera&);
	void setHomogeneousCorrespondences(EigMatX&, EigMatX&);

	// Get functions
	void getEssential(EigMat&);
	void getRotation(EigMat&);
	void getTranslation(EigVec&);
	void getEssentialResidual(scalar&);
	void getReconstructedPoints(EigMatX&);
    	void getCorrectedCorrespondences(EigMatX&, EigMatX&);
    	void getDataSize(int&);

	// Main functions
	void computeEssentialHL();
	void motionParametersKanatani();
    	void correctCorrespondencesKanatani();
	void triangulationKanatani();

	// Constructor / Destructor
	TwoImagesReconstruction();
	virtual ~TwoImagesReconstruction();
};

#endif // TWOIMAGESRECONSTRUCTION_H
