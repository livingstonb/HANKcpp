
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
	params.rho = 0.01;
	params.drs_N = 0.8;
	params.drs_Y = 0.9;
	params.setup(options);

	Model model = Model(params, income_dir);
	// std::cout << model.bgrid.size() << '\n';
	// std::cout << model.get_rb_effective().size() << '\n';

	SteadyState iss(model);

	HJB hjb(model, iss);
	// // hjb.maxiter = 1;
	hjb.iterate(iss);

	StationaryDist sdist;
	sdist.compute(model, iss, hjb);

	DistributionStatistics stats(params, model, hjb, sdist);

	std::cout << "E[NW] = " << stats.Enetworth << '\n';
}
