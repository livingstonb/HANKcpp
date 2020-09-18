
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
#include <ss_calibrator.h>

#include <minpack.h>

Parameters *global_params_ptr = NULL;
std::string income_dir = "2point_3_5";

SSCalibrator *global_calibrator_ptr = NULL;

void find_initial_steady_state(int n, double x[], double fvec[], int &iflag) {
	Parameters& p = *global_params_ptr;

	p.rho = exp(x[0]);
	p.rb = exp(x[p.nocc+2]);
	p.chi = x[p.nocc+3];
	p.update();

	Model model = Model(p, income_dir);

	SteadyState iss(p, model);
	iss.set(x, SteadyState::SSType::initial);
	iss.compute(SteadyState::SSType::initial);

	HJB hjb(model, iss);
	hjb.iterate(iss);

	StationaryDist sdist;
	sdist.gtol = 1.0e-9;
	sdist.compute(model, iss, hjb);

	DistributionStatistics stats(p, model, hjb, sdist);
	stats.print();

	// fvec[0] = stats.Ea / (iss.capital + iss.equity_A) - 1.0;
	// for (int io=0; io<p.nocc; ++io)
	// 	fvec[io+1] = stats.Elabor_occ[io] * model.occdist[io] / iss.labor_occ[io] - 1.0;

	// fvec[p.nocc+1] = stats.a_pctiles[5] / p.targetMedianIll - 1.0;
	// fvec[p.nocc+2] = stats.Eb / p.targetMeanLiq - 1.0;
	// fvec[p.nocc+3] = (stats.Ehours / p.hourtarget - 1.0) / 100.0;

	SSCalibrationArgs args(global_params_ptr, &model, &stats, &iss);
	global_calibrator_ptr->fill_fvec(args, fvec);
}

int main () {
	std::string income_dir = "2point_3_5";

	Options options; 
	options.fast = false;

	Parameters params;
	params.rho = 0.02;
	params.drs_N = 0.8;
	params.drs_Y = 0.9;
	params.dmax = 1e4;
	params.borrowing = true;
	// params.deathrate = 0.0;
	params.amax = 1000;
	params.na = 30;
	params.nb_pos = 30;
	// params.depreciation = 0.001;
	params.rb = 0.01 / 4.0;

	params.setup(options);
	global_params_ptr = &params;

	SSCalibrator calibrator(params);
	global_calibrator_ptr = &calibrator;

	Model model = Model(params, income_dir);

	// guess rho, chi,labor_occ, capital, and rb
	double x[params.nocc+4];
	x[0] = log(params.rho);
	for (int io=0; io<params.nocc; ++io)
		x[io+1] = log(params.hourtarget * params.meanlabeff * model.occdist[io]);

	x[params.nocc+1] = params.target_KY_ratio;
	x[params.nocc+2] = log(params.rb);
	x[params.nocc+3] = params.chi;

	int n = params.nocc + 4;
	double fvec[n];
	double tol = 1.0e-9;

	int lwa = n * (3 * n + 13);
	double wa[lwa];
	hybrd1(find_initial_steady_state, n, x, fvec, tol, wa, lwa);

	// int iflag=0;
	// find_initial_steady_state(n, x, fvec, iflag);
}