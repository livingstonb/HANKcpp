#ifndef _PROCEDURES_H
#define _PROCEDURES_H

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

std::vector<double> read_matrix(const std::string& file_loc);

std::size_t find_multiple(const std::string& line, int pos);

template<typename T>
void printvec(const T& vec) {
	for (auto x : vec) {
		std::cout << x << '\n';
	}
}

#endif