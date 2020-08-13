#ifndef _PROCEDURES_H
#define _PROCEDURES_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <mkl.h>
#include <mkl_cblas.h>
// #include <mkl_blas.h>
// #include <mkl_lapack.h>
// #include <mkl_lapacke.h>
#include <assert.h>
#include <hank.h>

template<typename T>
void linspace(double x, double y, int n, T& cont) {
	for (int i=0; i<n; ++i) {
		cont[i] = x + i * (y - x) / (n - 1);
	}
}

template<typename T>
void powerSpacedGrid(
	int n, double low, double high, double curv, T& grid)
{
	linspace(0.0, 1.0, n, grid);

	for (int i=0; i<n; ++i) {
		grid[i] = low + (high - low) * pow(grid[i], 1.0 / curv);
	}
}

template <typename T, typename F>
void apply(T& vec, F func) {
	std::for_each(vec.begin(), vec.end(), func);
}

template<typename T>
void printvec(const T& vec) {
	for (auto x : vec) {
		std::cout << x << '\n';
	}
}

template<typename T>
double vdot(const T& vec1, const T& vec2) {
	assert(vec1.size() == vec2.size());
	return cblas_ddot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
}

template<typename V>
std::pair<std::vector<double>,std::vector<double>> occupationGrid(const V& p)
{
	std::vector<double> occgrid, occdist;
	if (p.nocc == 1) {
		occgrid.push_back(0.0);
		occdist.push_back(1.0);
	}
	else {
		double lshareNY = (1.0 - p.alpha_N) * p.drs_N
			/ ((p.elast - 1.0) * (1.0 - p.alpha_Y) * p.drs_Y + (1.0 - p.alpha_N) * p.drs_N);

		if (lshareNY == 0.0) {
			occgrid.push_back(0.0);
			occdist.push_back(1.0 / p.nocc);
		}
		else if (lshareNY == 1.0) {
			occgrid.push_back(1.0);
			occdist.push_back(1.0 / p.nocc);
		}
		else {
			// Equally spaced in [0, 1], midpoints of intervals
			double lwidth = 1.0 / p.nocc;
			occgrid.resize(p.nocc);
			occgrid[0] = 0.5 * lwidth;
			occgrid[p.nocc-1] = 1.0 - 0.5 * lwidth;
			if (p.nocc >= 2) {
				for (int i=1; i<p.nocc-1; ++i) {
					occgrid[i] = occgrid[i-1] + lwidth;
				}
			}

			// Distribution has CDF x ^ par. Choose par to target mean a
			double lWNtoWY = 1.5;
			double lmeanocc = lshareNY / (lshareNY + (1.0 - lshareNY) * lWNtoWY);
			double lpar = lmeanocc / (1.0 - lmeanocc);
			for (int i=0; i<p.nocc; ++i) {
				occdist[i] = pow((occgrid[i] + 0.5 * lwidth), lpar)
					- pow((occgrid[i] - 0.5 * lwidth), lpar);
			}
		}
	}

	auto pair = std::make_pair(occgrid, occdist);
	return pair;
}

template<typename T>
double_vector vector2eigenv(const T& vec)
{
	double_vector out(vec.size());

	for (int i=0; i<vec.size(); ++i)
		out[i] = vec[i];	

	return out;
}

template<typename T>
double_matrix vector2eigenm(const T& vec, int n, int m)
{
	double_matrix out(n, m);

	assert(vec.size() == m * n);

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

template<typename T>
map_type boost2eigen(T& arr)
{
	int n0 = 1;
	int n1 = 1;
	int dims_found = 0;
	auto shape = arr.shape();

	for (int i=0; i<arr.num_dimensions(); ++i) {
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