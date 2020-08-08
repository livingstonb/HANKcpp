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
T linspace(double x, double y, int n) {
	T vec(n);

	for (int i=0; i<n; ++i) {
		vec[i] = x + i * (y - x) / (n - 1);
	}

	return vec;
}

template<typename T>
T powerSpacedGrid(
	int n, double low, double high, double curv)
{
	auto vec = linspace<T>(0.0, 1.0, n);

	for (auto it = vec.begin(); it != vec.end(); ++it) {
		*it = low + (high - low) * pow(*it, 1 / curv);
	}
	return vec;
}

template <typename T, typename F>
void apply(T& vec, F func) {
	std::for_each(vec.begin(), vec.end(), func);
}

template<typename T>
void printvec(const T& vec) {
	for (auto it = vec.begin(); it != vec.end(); ++it) {
		std::cout << *it << '\n';
	}
}

template<typename T>
double vdot(const T& vec1, const T& vec2) {
	assert(vec1.size() == vec2.size());
	return cblas_ddot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
}

template<typename T, typename V>
std::pair<T,T> occupationGrid(const V& p)
{
	T occgrid;
	T occdist;
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

template<size_t N>
array_type<N> new_array(const array_shape<N>& shape)
{
	array_type<N> arr(shape);
	// array_type<N> arr(shape, boost::fortran_storage_order());

	return arr;
}

template<typename T>
map_type to_eigen(T& arr)
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