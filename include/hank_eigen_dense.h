#ifndef HANK_EIGEN_DENSE
#define HANK_EIGEN_DENSE

#include <hank_config.h>
#include <Eigen/Core>
#include <hank_types.h>
#include <algorithm>

#ifdef HANK_INCLUDE_EIGEN_LU
#include <Eigen/LU>
#endif

using Eigen::seq;

using MatrixXr = Eigen::Matrix<hank_float_type, Eigen::Dynamic, Eigen::Dynamic>;

using VectorXr = Eigen::Matrix<hank_float_type, Eigen::Dynamic, 1>;

using ArrayXXr = Eigen::Array<hank_float_type, Eigen::Dynamic, Eigen::Dynamic>;

using ArrayXr = Eigen::Array<hank_float_type, Eigen::Dynamic, 1>;

using MapVectorXr = Eigen::Map<VectorXr>;

using MapMatrixXr = Eigen::Map<MatrixXr>;

using Eigen::VectorXi;

struct ModelMatrices {
	MatrixXr prodmarkov, ymarkov, ymarkovdiag, ymarkovoff;
};

inline VectorXr eflatten(const MatrixXr& arr) {
	const Eigen::Map<const VectorXr> arr_map(arr.data(), arr.size());
	VectorXr flattened = arr_map;
	return flattened;
}

inline VectorXr to_eigenv(const vector3dr& vec) {
	std::vector<hank_float_type> vcopy = vec.vector;
	VectorXr out = MapVectorXr(vcopy.data(), vcopy.size());
	return out;
}

template<typename T>
VectorXr vector2eigenv(const T& vec)
{
	VectorXr out(vec.size());

	for (size_t i=0; i<vec.size(); ++i)
		out[i] = vec[i];	

	return out;
}

template<typename T>
MatrixXr vector2eigenm(const T& vec, int n, int m)
{
	MatrixXr out(n, m);

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

template<typename V>
V as_eigen(const std::vector<double>& arr) {
	V out(arr.size());

	for (unsigned int i=0; i<arr.size(); ++i)
		out(i) = arr[i];

	return out;
}

template<typename V>
V as_eigen(const std::vector<long double>& arr) {
	V out(arr.size());

	for (unsigned int i=0; i<arr.size(); ++i)
		out(i) = arr[i];

	return out;
}

template<typename V>
V as_eigen(const vector3dr& arr) {
	return as_eigen<V>(arr.vector);
}

template<typename V>
V as_eigen(int n, const hank_float_type *arr) {
	V out(n);

	for (int i=0; i<n; ++i)
		out(i) = arr[i];

	return out;
}

template<typename T, typename V>
Eigen::Map<T> as_eigen_map(V& arr) {
	return Eigen::Map<T>(arr.data(), arr.size());
}

template<typename T>
VectorXr cumsum(const T& arr) {
	VectorXr out(arr.size());
	double val = 0.0;

	for (unsigned int i=0; i<arr.size(); ++i) {
		val += arr(i);
		out[i] = val;
	}

	return out;
}

inline MatrixXr deye(int n) {
	MatrixXr mat = MatrixXr::Zero(n, n);

	for (int i=0; i<n; ++i)
		mat(i, i) = 1;

	return mat;
}

namespace EigenFunctions {
	template<typename T>
	double sum(const T& arr) {
		auto earr = as_eigen_map<const VectorXr>(arr);
		return earr.sum();
	}

	template<typename T, typename V>
	double dot(const T& arr1, const V& arr2) {
		auto earr1 = as_eigen_map<const VectorXr>(arr1);
		auto earr2 = as_eigen_map<const VectorXr>(arr2);

		return earr1.dot(earr2);
	}
}


#endif