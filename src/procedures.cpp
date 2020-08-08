#include <procedures.h>

std::vector<double> linspace(double x, double y, int n) {
	std::vector<double> vec(n);

	for (int i=0; i<n; ++i) {
		vec[i] = x + i * (y - x) / (n - 1);
	}

	return vec;
}

std::vector<double> PowerSpacedGrid(
	int n, double low, double high, double curv)
{
	auto vec = linspace(0.0, 1.0, n);

	for (auto it = vec.begin(); it != vec.end(); ++it) {
		*it = low + (high - low) * pow(*it, 1 / curv);
	}
	return vec;
}

