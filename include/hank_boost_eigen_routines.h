#ifndef _HANK_BOOST_EIGEN_ROUTINES_H
#define _HANK_BOOST_EIGEN_ROUTINES_H

#include <hank_config.h>
#include <hank_eigen_dense.h>

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
map_type boost2eigen(T& arr)
{
	int n0 = 1;
	int n1 = 1;
	int dims_found = 0;
	auto shape = arr.shape();

	for (size_t i=0; i<arr.num_dimensions(); ++i) {
		if (shape[i] > 1) {
			if (dims_found == 0) {
				n0 = shape[i];
			}
			else if (dims_found == 1) {
				n1 = shape[i];
			}

			++dims_found;
		}
	}

	// Cannot squeeze array into two dimensions
	assert(dims_found <= 2);

	map_type map(arr.data(), n0, n1);
	return map;
}

template<typename T, typename V>
double boost_dot(const T& boost_arr, const V& eigen_arr) {
	// auto flat = flatten_array3d(boost_arr);
	// double dotprod = 0.0;

	// for (int i=0; i<flat.num_elements(); ++i)
	// 	dotprod += flat[i][0][0] * eigen_arr[i];
	double dotprod = 0.0;

	auto shape = boost_arr.shape;
	for (int i0=0; i0<shape[0]; ++i0)
		for (int i1=0; i1<shape[1]; ++i1)
			for (int i2=0; i2<shape[2]; ++i2)
				dotprod += boost_arr(i0, i1, i2) * eigen_arr[i0 + shape[0] * i1 + shape[0] * shape[1] * i2];

	return dotprod;
}

#endif