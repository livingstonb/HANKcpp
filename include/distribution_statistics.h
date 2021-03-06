#ifndef _DISTRIBUTION_STATISTICS_H
#define _DISTRIBUTION_STATISTICS_H

#include <hank_config.h>
#include <hank.h>
#include <vector>
#include <map>

// Forward declarations
class Parameters;

class Model;

class StationaryDist;

namespace Upwinding {
	class Policies;
}

// Class for distribution statistics
class DistributionStatistics : public HankBase
{
	public:
		DistributionStatistics() {}

		DistributionStatistics(const Parameters& p, const Model& model,
			const Upwinding::Policies& policies, const vector3dr& density);

		DistributionStatistics& operator=(const DistributionStatistics* other_stats)
		{
			*this = other_stats;
			return *this;
		}

		std::string title() const override {return "STATISTICS";}

		std::map<std::string, hank_float_type> variables_map() const override;

		vector3dr density;

		double Ehours, Ea, Eb, Enetworth, Elabor;

		std::vector<double> Elabor_occ, a_pctiles, b_pctiles, nw_pctiles;

		std::vector<double> pctiles = {0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 0.99};
};

#endif