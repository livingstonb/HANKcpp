#ifndef _DISTRIBUTION_STATISTICS_H
#define _DISTRIBUTION_STATISTICS_H

// Forward declarations
class Parameters;

class Model;

class HJB;

class StationaryDist;

namespace Upwinding {
	class Policies;
}

template<typename T, typename V>
double boost_dot(const T& boost_arr, const V& eigen_arr) {
	double dotprod = 0.0;

	auto shape = boost_arr.shape;
	for (int i0=0; i0<shape[0]; ++i0)
		for (int i1=0; i1<shape[1]; ++i1)
			for (int i2=0; i2<shape[2]; ++i2)
				dotprod += boost_arr(i0, i1, i2) * eigen_arr[i0 + shape[0] * i1 + shape[0] * shape[1] * i2];

	return dotprod;
}

// Class for distribution statistics
class DistributionStatistics {

	class DistStruct;

	public:
		DistributionStatistics(const Parameters& p, const Model& model,
			const HJB& hjb, const StationaryDist& sdist);

		void compute_moments(const Upwinding::Policies& policies, const DistStruct& grids);

		void print();

		double Ehours, Ea, Eb, Enetworth;
};

#endif