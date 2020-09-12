#ifndef _DISTRIBUTION_STATISTICS_H
#define _DISTRIBUTION_STATISTICS_H

// Forward declarations
class Parameters;

class Model;

class HJB;

class StationaryDist;

class DistGrids;

namespace Upwinding {
	class Policies;
}

// Class for distribution statistics
class DistributionStatistics {
	public:
		DistributionStatistics(const Parameters& p, const Model& model,
			const HJB& hjb, const StationaryDist& sdist);

		void compute_moments(const Model& model,
			const Upwinding::Policies& policies, const DistGrids& grids);

		double Ehours, Ea, Eb, Enetworth;
};

#endif