#ifndef HANK_EIGEN_DENSE
#define HANK_EIGEN_DENSE

#include <hank_config.h>
#include <Eigen/Core>
#include <hank_types.h>

using Eigen::seq;

typedef Eigen::Map<Eigen::MatrixXd> map_type;

typedef Eigen::Map<Eigen::VectorXd> map_type_vec;

typedef Eigen::VectorXd double_vector;

typedef Eigen::ArrayXd double_array;

typedef Eigen::MatrixXd double_matrix;

typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> bool_vector;

typedef double_vector grid_type;

inline double_vector eflatten(const double_matrix& arr) {
	const Eigen::Map<const double_vector> arr_map(arr.data(), arr.size());
	double_vector flattened = arr_map;
	return flattened;
}

inline double_vector to_eigenv(const std::vector<double>& vec) {
	std::vector<double> vcopy = vec;
	double_vector out = map_type_vec(vcopy.data(), vcopy.size());
	return out;
}

inline double_vector to_eigenv(const StdVector3d<double>& vec) {
	std::vector<double> vcopy = vec.vector;
	double_vector out = map_type_vec(vcopy.data(), vcopy.size());
	return out;
}

#endif