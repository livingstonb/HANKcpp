
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
	SSCalibrator& cal = *global_calibrator_ptr;

	std::cout << "\nCalibration parameters updated:\n";

	cal.update_params(&p, x);

	Model model = Model(p, income_dir);

	SteadyState iss(p, model);
	cal.update_ss(p, &iss, x);
	std::cout << '\n';

	iss.compute(SteadyState::SSType::initial);

	HJB hjb(model, iss);
	hjb.iterate(iss);

	StationaryDist sdist;
	sdist.gtol = 1.0e-9;
	sdist.compute(model, iss, hjb);

	DistributionStatistics stats(p, model, hjb, sdist);
	stats.print();

	SSCalibrationArgs args(global_params_ptr, &model, &stats, &iss);
	global_calibrator_ptr->fill_fvec(args, fvec);

	VectorXd fvals(n);
	for (int i=0; i<n; ++i)
		fvals[i] = fvec[i];

	std::cout << "fvec = \n";
	for (int i=0; i<n; ++i)
		std::cout << "  " << fvec[i] << '\n';
	std::cout << "  fvec norm = " << fvals.norm() << '\n';
	std::cout << "--------------------------\n\n";
}

int main () {
	std::string income_dir = "2point_3_5";

	Options options; 
	options.fast = false;

	Parameters params;
	params.rho = 0.015;
	params.drs_N = 0;
	params.drs_Y = 0.9;
	params.dmax = 1e3;
	params.borrowing = true;
	// params.deathrate = 0.0;
	params.amax = 500;
	params.na = 50;
	params.nb_pos = 30;
	params.depreciation = 0.05 / 4;
	params.elast = 2;
	params.nocc = 1;
	// params.riskaver = 1.0;
	params.rb = 0.02 / 4.0;

	params.setup(options);
	global_params_ptr = &params;

	SSCalibrator calibrator(params);
	global_calibrator_ptr = &calibrator;

	Model model = Model(params, income_dir);

	// guess rho, chi,labor_occ, capital, and rb
	double x[params.nocc+4];
	calibrator.fill_xguess(params, model, x);

	int n = params.nocc + 4;
	double fvec[n];
	double tol = 1.0e-9;

	assert( n == calibrator.nmoments() );

	int lwa = n * (3 * n + 13);
	double wa[lwa];
	hybrd1(find_initial_steady_state, n, x, fvec, tol, wa, lwa);

	// int iflag=0;
	// find_initial_steady_state(n, x, fvec, iflag);
}