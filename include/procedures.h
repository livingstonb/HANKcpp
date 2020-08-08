#ifndef _PROCEDURES_H
#define _PROCEDURES_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

std::vector<double> linspace(double x, double y, int n);

std::vector<double> PowerSpacedGrid(
	int n, double low, double high, double curv);

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

#endif