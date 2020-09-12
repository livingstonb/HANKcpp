#ifndef _DISTRIBUTION_STATISTICS_H
#define _DISTRIBUTION_STATISTICS_H

#include <vector>
#include <hank_eigen_dense.h>

class Parameters;

class Model;

class HJB;

class StationaryDist;

class DistGrids {
	public:
		double_matrix networth;
};

namespace Upwinding {
	class Policies;
}

class DistributionStatistics {
	public:
		DistributionStatistics(const Parameters& p, const Model& model,
			const HJB& hjb, const StationaryDist& sdist);

		void compute_moments(const Model& model,
			const Upwinding::Policies& policies, const DistGrids& grids);

		DistGrids grids;
		double_matrix pmass;
		double_vector pmass1d;
		double Ehours;
		double Ea, Eb, Enetworth;
};

#endif