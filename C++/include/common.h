//
// Created by Dimitris Zermas on 4/8/17.
//

#ifndef __COMMON_H_INCLUDED__
#define __COMMON_H_INCLUDED__

#define PI 3.14

// You can change the types of the Eigen matrices and vectors here
// Most common change would be for scalar to become double from float
// This header file is seen by all other files
typedef float scalar;
typedef Eigen::Matrix<scalar,3,3> EigMat;
typedef Eigen::Matrix<scalar,3,1> EigVec;
typedef Eigen::Matrix<scalar,Eigen::Dynamic,Eigen::Dynamic> EigMatX;
typedef Eigen::Matrix<scalar,Eigen::Dynamic,1> EigVecX;

/* Generates random numbers in the [0,1] range */
template <typename T>
inline T rand01(void){
    return std::rand()/(T)RAND_MAX;
}

/* Generates a "random" rotation. The rotation matrix is constrained so that it makes sense for this application */
template <typename Mat, typename Scalar>
inline Mat generateRandomRotation(void){

    Scalar rand_angle_1 = rand01<scalar>() * 5.0 * PI / 180.0;
    Scalar rand_angle_2 = rand01<scalar>() * -5.0 * PI / 180.0;
    Scalar rand_angle_3 = rand01<scalar>() * 5.0 * PI / 180.0;

    Mat Rz;
	Rz << std::cos(rand_angle_3), -std::sin(rand_angle_3), 0,
		  std::sin(rand_angle_3),  std::cos(rand_angle_3), 0,
		  0						,  0					 , 1;

    Mat Ry;
	Ry <<  std::cos(rand_angle_2), 0, std::sin(rand_angle_2),
		   0					 , 1, 0,
		  -std::sin(rand_angle_2), 0, std::cos(rand_angle_2);


    Mat Rx;
	Rx << 1, 0					   ,  0,	 
		  0, std::cos(rand_angle_1), -std::sin(rand_angle_1),
		  0, std::sin(rand_angle_1),  std::cos(rand_angle_1);

    return Rz * Ry * Rx;
}

/* Skew Symmetric matrix from a vector */
template <typename Mat, typename Vec>
inline Mat skewSymmetric(Vec &v){
    Mat V;
	V << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;
	return V;
}

/* Convert std::vector 3D data to Eigen format */
template <typename T, typename Mat>
inline void vec2eigen(std::vector<T> &data, Mat &C1_P){
    int n = data.size() / 3;
    C1_P.resize(n,3);
    for (size_t i = 0; i < data.size(); i++){
        C1_P(i/3,i%3) = data.at(i);
    }
}

/* Read the .txt file that contains the 3D points */
template <typename Mat>
inline void readData(const std::string &filename, Mat &data_eigen){
    std::ifstream in;
    in.open(filename);
    if (!in){
        std::cout << "Error opening file" << std::endl;
        exit(1);
    }

    scalar x;
    std::vector<scalar> data;
    while (in >> x){
        data.push_back(x);
    }

    in.close();

    vec2eigen(data, data_eigen); // Convert data to Eigen format
}

/* Generic print for Eigen Matrices. Eigen makes sure the matrix is valid */
template <typename Mat>
inline void printMatrix(Mat& V){
    int n = V.rows();
    int m = V.cols();
	for (size_t i=0;i<n;i++){
		for (size_t j=0;j<m;j++){
			std::cout << V(i,j) << "\t";
		}
		std::cout << std::endl;
	}
}

/* Generic print for Eigen Vectors. Eigen makes sure the vector is valid */
template <typename Vec>
inline void printVector(Eigen::MatrixBase<Vec>& V){
    int n = V.rows();
    int m = V.cols();

    if (n != 1 && m != 1){
        std::cout << "printVector: input is not a vector" << std::endl;
        exit(1);
    }
    // vertical or horizontal vector?
    if (n == 1){
        for (size_t i=0;i<m;i++){
            std::cout << V(i) << "\t";
        }
    }
    else{
        for (size_t i=0;i<n;i++){
            std::cout << V(i) << "\n";
        }
    }
    std::cout << std::endl;
}

/* Computes the norm L2 of an Eigen vector */
template <typename Scalar, typename Vec>
inline Scalar norm(Vec& V){
    int n = V.rows();
    int m = V.cols();

    if (n != 1 && m != 1){
        std::cout << "norm: input is not a vector" << std::endl;
        exit(1);
    }

    Scalar d = 0.0;
    int v_size = (n > m) ? n : m;
    for (size_t i=0;i<v_size;i++){
        d += V(i) * V(i);
    }
    d = std::sqrt(d);

    return d;
}

/* Computes the dot product of two Eigen-type vectors */
template <typename Scalar, typename Vec>
inline Scalar dot(Vec& V, Vec& U){
    int n = V.rows();
    int m = V.cols();
    int k = U.rows();
    int l = U.cols();

    if (n != 1 && m != 1 && k != 1 && l != 1){
        std::cout << "dot: at least one input is not a vector" << std::endl;
        exit(1);
    }

    int v_size = (n > m) ? n : m;
    int u_size = (k > l) ? k : l;
    if (v_size != u_size){
        std::cout << "dot: two inputs not of the same size" << std::endl;
        exit(1);
    }

    Scalar d = 0.0;
    for (size_t i=0;i<v_size;i++){
        d += V(i) * U(i);
    }

    return d;
}

/* Computes the cross product of two Eigen-type vectors */
template <typename Vec>
inline Vec cross(Vec& V, Vec& U){
    int n = V.rows();
    int m = V.cols();
    int k = U.rows();
    int l = U.cols();

    if (n != 1 && m != 1 && k != 1 && l != 1){
        std::cout << "cross: at least one input is not a vector" << std::endl;
        exit(1);
    }

    int v_size = (n > m) ? n : m;
    int u_size = (k > l) ? k : l;
    if (v_size != u_size){
        std::cout << "cross: two inputs not of the same size" << std::endl;
        exit(1);
    }

    if (v_size != 3){
        std::cout << "cross: inputs not of size 3" << std::endl;
        exit(1);
    }

    Vec vxu;
    vxu(0) = V(1) * U(2) - U(1) * V(2);
    vxu(1) = U(0) * V(2) - V(0) * U(2);
    vxu(2) = V(0) * U(1) - U(0) * V(1);

    return vxu;
}

/* Triple product as defined by Kenichi Kanatani's books:
   (scalar) triple(x1,x2,x3) = dot(cross(x1,x2),x3) */
template <typename Scalar, typename Vec>
inline Scalar triple(Vec& V, Vec& U, Vec& W){
    Scalar product;
    Vec c = cross<Vec>(V,U);
    product = dot<Scalar,Vec>(c,W);
    return product;
}

/* Computes the normalized vector of an Eigen-type vector */
template <typename Vec>
inline Vec normalize(Vec& V){
    scalar n = norm<scalar,Vec>(V);
    Vec normalized;
    normalized(0) = V(0) / n;
    normalized(1) = V(1) / n;
    normalized(2) = V(2) / n;
    return normalized;
}

/* Returns the sign of a value */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/* Solve homogeneous equation */
template <typename Mat, typename Vec>
inline void solveHomogeneousEq(Mat& A, Vec& x){
    Eigen::JacobiSVD<EigMatX> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    EigMatX V = svd.matrixV();
    x = V.col(V.cols()-1); // last column of V is the solution
};

#endif //COMMON_H
