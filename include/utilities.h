#ifndef _PROCEDURES_H
#define _PROCEDURES_H

#include <hank_config.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <assert.h>
#include <fstream>

namespace HankUtilities {

std::vector<double> read_matrix(const std::string& file_loc);

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

void check_cminpack_success(int info);

}

#endif