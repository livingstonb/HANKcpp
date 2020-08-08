#ifndef _PROCEDURES_H
#define _PROCEDURES_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#include <mkl_lapacke.h>
#include <assert.h>

std::vector<double> linspace(double x, double y, int n);

std::vector<double> powerSpacedGrid(
	int n, double low, double high, double curv);

// template<typename T>
// auto occupationGrid(const T& p);

template <typename F>
void apply(std::vector<double>& vec, F func) {
	std::for_each(vec.begin(), vec.end(), func);
}

template<typename T>
void printvec(T& vec) {
	for (auto it = vec.begin(); it != vec.end(); ++it) {
		std::cout << *it << '\n';
	}
}

template<typename T>
double vdot(T vec1, T vec2) {
	assert(vec1.size() == vec2.size());
	return cblas_ddot(vec1.size(), vec1.data(), 1, vec2.data(), 1);
}

template<typename T>
auto occupationGrid(const T& p)
{
	std::vector<double> occgrid;
	std::vector<double> occdist;
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

#endif