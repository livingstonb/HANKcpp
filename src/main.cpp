
#include <iostream>
#include <parameters.h>
#include <hank_types.h>
#include <distribution_statistics.h>
#include <model.h>
#include <steady_state.h>
#include <bellman.h>
#include <stationary_dist.h>
#include <adjustment_costs.h>
#include <utilities.h>

int main () {
	std::string income_dir = "2point_3_5";

	Options options;
	options.fast = false;

	Parameters params;
	params.rho = 0.02;
	params.drs_N = 0.8;
	params.drs_Y = 0.9;
	// params.na = 25;
	// params.nb_pos = 25;
	// params.depreciation = 0.001;
	params.rb = 0.001 / 4.0;
	params.setup(options);

	Model model = Model(params, income_dir);

	SteadyState iss(model);

	HJB hjb(model, iss);
	hjb.iterate(iss);

	StationaryDist sdist;
	sdist.compute(model, iss, hjb);

	DistributionStatistics stats(params, model, hjb, sdist);
	stats.print();
}
