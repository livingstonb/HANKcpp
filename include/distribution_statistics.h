#ifndef _DISTRIBUTION_STATISTICS_H
#define _DISTRIBUTION_STATISTICS_H

#include <vector>

// Forward declarations
class Parameters;

class Model;

class HJB;

class StationaryDist;

namespace Upwinding {
	class Policies;
}

// Class for distribution statistics
class DistributionStatistics {

	class DistStruct;

	public:
		DistributionStatistics() {}

		DistributionStatistics(const Parameters& p, const Model& model,
			const HJB& hjb, const StationaryDist& sdist);

		DistributionStatistics& operator=(const DistributionStatistics* other_stats) {
			*this = other_stats;
			return *this;
		}

		void print();

		double Ehours, Ea, Eb, Enetworth, Elabor;

		std::vector<double> Elabor_occ, a_pctiles, b_pctiles, nw_pctiles;

		std::vector<double> pctiles = {0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 0.99};
};

#endif