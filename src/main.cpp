
#include <iostream>
#include <hank_config.h>
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

const Options *global_hank_options = NULL;

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

	cal.print_fvec(fvec);
}

void set_to_fortran_params(Parameters& p) {

}

int main () {
	std::string income_dir = "2point_3_5";

	Options options; 
	options.fast = false;
	options.print_diagnostics = true;

	global_hank_options = &options;

	// Parameters
	Parameters params;
	params.rho = 0.022;
	params.drs_N = 0;
	params.drs_Y = 0.9;
	params.dmax = 1e6;
	params.borrowing = true;
	// params.deathrate = 0.0;
	params.amax = 400;
	params.na = 40;
	params.nb_pos = 40;
	params.depreciation = 0.05 / 4;
	params.elast = 2;
	params.nocc = 1;
	// params.riskaver = 1.0;
	params.rb = 0.03 / 4.0;

	params.setup(options);
	global_params_ptr = &params;

	// Calibration
	SSCalibrator calibrator;
	calibrator.calibrateLaborDisutility = false;
	calibrator.calibrateRb = true;
	calibrator.calibrateDiscountRate = true;
	calibrator.setup(params);
	global_calibrator_ptr = &calibrator;

	Model model = Model(params, income_dir);

	// guess rho, chi,labor_occ, capital, and rb
	int n = calibrator.nmoments;
	double x[n];
	calibrator.fill_xguess(params, model, x);

	double fvec[n];
	double tol = 1.0e-9;

	int lwa = n * (3 * n + 13);
	double wa[lwa];
	hybrd1(find_initial_steady_state, n, x, fvec, tol, wa, lwa);

	// int iflag=0;
	// find_initial_steady_state(n, x, fvec, iflag);
}