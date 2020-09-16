
#include <iostream>
#include <parameters.h>
#include <hank_types.h>
#include <model.h>
#include <steady_state.h>
#include <bellman.h>
#include <stationary_dist.h>
#include <distribution_statistics.h>
#include <adjustment_costs.h>
#include <utilities.h>
#include <math.h>

#include <minpack.h>

Parameters *global_params_ptr = NULL;

Model *global_model_ptr = NULL;

SteadyState *global_iss_ptr = NULL;

// void find_initial_steady_state(int n, double *x, double *fvec, int *iflag) {
// 	SteadyState& iss = *global_iss_ptr;

// 	global_iss_ptr->set(x, SteadyState::SSType::initial);
// 	global_iss_ptr->compute(SteadyState::SSType::initial);

// 	HJB hjb(*global_model_ptr, *global_iss_ptr);
// 	hjb.iterate(*global_iss_ptr);

// 	StationaryDist sdist;
// 	sdist.compute(*global_model_ptr, *global_iss_ptr, hjb);

// 	DistributionStatistics stats(*global_params_ptr, *global_model_ptr, hjb, sdist);
// 	stats.print();

// 	fvec[0] = stats.Ea / (iss.capital + iss.equity_A) - 1.0;
// 	for (int i=1; i<global_params_ptr->nocc+4; ++i)
// 		fvec[i] = 0.2;
// }

void find_initial_steady_state(int n, double x[], double fvec[], int &iflag) {
	SteadyState& iss = *global_iss_ptr;

	global_iss_ptr->set(x, SteadyState::SSType::initial);
	global_iss_ptr->compute(SteadyState::SSType::initial);

	HJB hjb(*global_model_ptr, *global_iss_ptr);
	hjb.iterate(*global_iss_ptr);

	StationaryDist sdist;
	sdist.compute(*global_model_ptr, *global_iss_ptr, hjb);

	DistributionStatistics stats(*global_params_ptr, *global_model_ptr, hjb, sdist);
	stats.print();

	fvec[0] = stats.Ea / (iss.capital + iss.equity_A) - 1.0;
	for (int i=1; i<global_params_ptr->nocc+4; ++i)
		fvec[i] = 0.2;
}

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
	global_params_ptr = &params;

	Model model = Model(params, income_dir);
	global_model_ptr = &model;

	// guess rho, chi,labor_occ, capital, and rb
	double chi = 0.5;
	double x[params.nocc+4];
	x[0] = log(params.rho);
	for (int io=0; io<params.nocc; ++io)
		x[io+1] = params.hourtarget * params.meanlabeff * model.occdist[io];

	x[params.nocc+1] = params.target_KY_ratio;
	x[params.nocc+2] = log(params.rb);
	x[params.nocc+3] = chi;

	SteadyState iss(params, model);
	global_iss_ptr = &iss;

	// iss.set(x, SteadyState::SSType::initial);
	// iss.compute(SteadyState::SSType::initial);

	
	// HJB hjb(model, iss);
	// hjb.iterate(iss);

	// StationaryDist sdist;
	// sdist.compute(model, iss, hjb);

	// DistributionStatistics stats(params, model, hjb, sdist);
	// stats.print();

	int n = params.nocc + 4;
	double fvec[n];
	double tol = 1.0e-9;

	int lwa = n * (3 * n + 13);
	double wa[lwa];

	int iflag = 0;
	// find_initial_steady_state(1, x, fvec, iflag);

	hybrd1(find_initial_steady_state, n, x, fvec, tol, wa, lwa);
}