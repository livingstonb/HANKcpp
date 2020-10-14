#ifndef HANK_EIGEN_SPARSE
#define HANK_EIGEN_SPARSE

#include <hank_config.h>
#include <Eigen/SparseCore>
#include <vector>

using SparseXd = Eigen::SparseMatrix<double>;

struct SparseMatContainer {
	SparseMatContainer(SparseXd matrix_) : matrix(matrix_) {}

	SparseXd matrix;

	SparseXd& get() {return matrix;}
};

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

using triplet_type = Eigen::Triplet<hank_float_type>;

using triplet_list = std::vector<triplet_type>;

inline SparseXd speye(int n) {
	SparseXd mat(n, n);
	triplet_list trips;
	trips.reserve(n);

	for (int i=0; i<n; ++i)
		trips.push_back(triplet_type(i, i, 1.0));

	mat.setFromTriplets(trips.begin(), trips.end());
	return mat;
}

#endif