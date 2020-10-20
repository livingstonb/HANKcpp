#ifndef _PROCEDURES_H
#define _PROCEDURES_H

#include <hank_config.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <fstream>

namespace HankUtilities {

std::vector<hank_float_type> read_matrix(const std::string& file_loc);

std::size_t find_multiple(const std::string& line, int pos);

template<typename T>
void printvec(const T& vec) {
	for (auto x : vec) {
		std::cout << x << '\n';
	}
}

template<typename T>
void printvec(const T* vec, int n) {
	for (int i=0; i<n; ++i) {
		std::cout << vec[i] << '\n';
	}
}

inline void horzline() {
	std::cout << "----------------------------\n";
}

inline void print_values(const std::vector<std::string>& names, const std::vector<double>& values) {
	for (unsigned int i=0; i<names.size(); ++i) {
		std::cout << "  " << names[i] << " = " << values[i] << '\n';
	}
}

template<typename T, typename V>
void fillarr(T* arr, V val, int n) {
	for (int i=0; i<n; ++i)
		arr[i] = val;
}

template<typename T, typename V>
void fillarr(T* arr, V val, int n, int m) {
	for (int i=0; i<n; ++i)
		for (int j=0; j<m; ++j)
			arr[i + n * j] = val;
}

template<typename T, typename V>
void copyarr(const T* arr_in, V* arr_out, int n) {
	for (int i=0; i<n; ++i)
		arr_out[i] = arr_in[i];
}

}

#endif