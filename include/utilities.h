#ifndef _PROCEDURES_H
#define _PROCEDURES_H

#include <hank.h>

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
// #include <mkl.h>
// #include <mkl_cblas.h>
// #include <mkl_blas.h>
// #include <mkl_lapack.h>
// #include <mkl_lapacke.h>
#include <assert.h>
#include <fstream>

class Parameters;

template<typename T>
void linspace(double x, double y, int n, T& cont) {
	for (int i=0; i<n; ++i) {
		cont[i] = x + i * (y - x) / (n - 1);
	}
}

void powerSpacedGrid(double low, double high, double curv, grid_type& grid);

void adjustPowerSpacedGrid(grid_type& grid);

std::vector<double> read_matrix(const std::string& file_loc);

std::size_t find_multiple(const std::string& line, int pos);

sparse_matrix speye(int n);

template<typename T>
void printvec(const T& vec) {
	for (auto x : vec) {
		std::cout << x << '\n';
	}
}

// template<typename T>
// double vdot(const T& vec1, const T& vec2) {
// 	assert(vec1.size() == vec2.size());
// 	return cblas_ddot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
// }


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


template<typename T, size_t N>
boost_array_type<T, N> new_array(const boost_array_shape<T, N>& shape)
{
	boost_array_type<T, N> arr(shape);
	// array_type<N> arr(shape, boost::fortran_storage_order());

	return arr;
}

template<typename T, size_t N>
boost_array_type<T, N> reshape_array(const boost_array_type<T, N>& arr, const boost_array_shape<T, N>& shape)
{
	boost_array_type<T, N> out = arr;
	out.reshape(shape);
	// boost_array_shape<T, N> new_shape;
	// int n = arr.num_dimensions();

	// if (n > 1) {
	// 	auto new_shape = arr.shape();
	// 	new_shape[0] = arr.num_elements();

	// 	for (int i=1; i<n; ++i)
	// 		new_shape[i] = 1;

	// 	out.reshape(new_shape);
	// }

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

#endif