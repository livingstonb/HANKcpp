#ifndef HANK_EIGEN_DENSE
#define HANK_EIGEN_DENSE

#include <hank_config.h>
#include <Eigen/Core>
#include <hank_types.h>
#include <algorithm>

using Eigen::seq;

typedef Eigen::Map<Eigen::MatrixXd> map_type;

typedef Eigen::Map<Eigen::VectorXd> map_type_vec;

typedef Eigen::VectorXd double_vector;

typedef Eigen::ArrayXXd double_array;

typedef Eigen::ArrayXd double_array_vector;

typedef Eigen::MatrixXd double_matrix;

typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> bool_vector;

typedef double_vector grid_type;

using Eigen::ArrayXd;

using Eigen::VectorXd;

using Eigen::MatrixXd;

using Eigen::VectorXi;

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

template<typename T>
double_vector vector2eigenv(const T& vec)
{
	double_vector out(vec.size());

	for (size_t i=0; i<vec.size(); ++i)
		out[i] = vec[i];	

	return out;
}

template<typename T>
double_matrix vector2eigenm(const T& vec, int n, int m)
{
	double_matrix out(n, m);

	assert(vec.size() == static_cast<size_t>(m * n));

	for (int i=0; i<n; ++i)
		for (int j=0; j<m; ++j)
			out(i, j) = vec[i*m + j];	

	return out;
}

template<typename T>
std::vector<double> to_vector(const T& emat) {
	std::vector<double> vec(emat.size());
	std::copy(emat.begin(), emat.end(), vec.begin());
	return vec;
}

template<typename T>
map_type_vec to_eigen(T& arr) {
	map_type_vec map(arr.data(), arr.size());
	return map;
}

template<typename T>
map_type to_eigen(T& arr, int n, int m) {
	map_type map(arr.data(), n, m);
	return map;
}

template<typename V>
V as_eigen(std::vector<double>& arr) {
	V out(arr.size());

	for (unsigned int i=0; i<arr.size(); ++i)
		out(i) = arr[i];

	return out;
}

template<typename V>
V as_eigen(StdVector3d<double>& arr) {
	return as_eigen<V>(arr.vector);
}

#endif