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