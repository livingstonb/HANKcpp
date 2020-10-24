#ifndef HANK_EIGEN_SPARSE
#define HANK_EIGEN_SPARSE

#include <hank_config.h>
#include <Eigen/SparseCore>
#include <vector>

using SparseXd = Eigen::SparseMatrix<double>;

struct SparseMatContainer {
	SparseMatContainer(SparseXd&& matrix_) : matrix(matrix_) {}

	SparseXd matrix;

	SparseXd& get() {return matrix;}
};

// struct EigenTripletContainer {
// 	EigenTripletContainer()
// }

#if HANK_EIGEN_SPARSE_SOLVER == 0
	#include <Eigen/SparseQR>
	using sparse_solver = Eigen::SparseQR<SparseXd, Eigen::COLAMDOrdering<int>>;
#elif HANK_EIGEN_SPARSE_SOLVER == 1
	#include <Eigen/SparseLU>
	using sparse_solver = Eigen::SparseLU<SparseXd, Eigen::COLAMDOrdering<int>>;
#elif HANK_EIGEN_SPARSE_SOLVER == 2
	#include <Eigen/UmfPackSupport>
	using sparse_solver = Eigen::UmfPackLU<SparseXd>;
#endif

using EigenTriplet = Eigen::Triplet<hank_float_type>;

inline SparseXd speye(int n) {
	SparseXd mat(n, n);
	std::vector<EigenTriplet> trips;
	trips.reserve(n);

	for (int i=0; i<n; ++i)
		trips.push_back(EigenTriplet(i, i, 1.0));

	mat.setFromTriplets(trips.begin(), trips.end());
	return mat;
}

#endif