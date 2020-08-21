#ifndef HANK_EIGEN_SPARSE
#define HANK_EIGEN_SPARSE

#define NDEBUG
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <vector>

using sparse_matrix = Eigen::SparseMatrix<double>;

using sparse_solver = Eigen::SparseLU<sparse_matrix, Eigen::COLAMDOrdering<int>>;

using triplet_type = Eigen::Triplet<double>;

using triplet_list = std::vector<triplet_type>;

sparse_matrix speye(int n) {
	sparse_matrix mat(n, n);
	triplet_list trips;
	trips.reserve(n);

	for (int i=0; i<n; ++i)
		trips.push_back(triplet_type(i, i, 1.0));

	mat.setFromTriplets(trips.begin(), trips.end());
	return mat;
}

#endif