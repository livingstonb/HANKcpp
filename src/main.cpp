
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

#include <cminpack.h>
#include <cminpackP.h>

const Options *global_hank_options = NULL;

Parameters *global_params_ptr = NULL;
std::string income_dir = "2point_3_5";

SSCalibrator *global_calibrator_ptr = NULL;

int fcn(void* /* _p */, int n, const real *x, real *fvec, int /* iflag */ ) {
	Parameters& p = *global_params_ptr;
	SSCalibrator& cal = *global_calibrator_ptr;

	assert(cal.nmoments == n);

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

	return 0;
}

void set_to_fortran_params(Parameters& p) {
	p.na = 40;
	p.nb_pos = 40;
	p.nb_neg = 10;
	p.nocc = 1;
	p.bcurv = 0.35;
	p.bcurv_neg = 0.4;
	p.acurv = 0.15;
	p.amax = 2000;
	p.bmax = 40;
	p.cmin = 1e-5;
	p.dmax = 100;
	p.imposeMaxHours= true;
	p.perfectAnnuityMarkets = true;
	p.adjustProdGridFrisch = true;
	p.adjFrischGridFrac = 0.85;
	p.taxHHProfitIncome = true;
	p.profdistfracA = 1;
	p.profdistfracB = 0;
	p.profdistfracW = 0;
	p.profdistfracL = 0;
	p.rho = 0.055 / 4.0;
	p.deathrate = 1.0 / (4.0 * 45.0);
	p.riskaver = 1.5;
	p.rb = 0.015 / 4.0;
	p.borrwedge = 0.025;
	p.blim = -1.0;
	p.frisch = 1.0;
	p.hourtarget = 1.0 / 3.0;
	p.drs_Y = 1.0;
	p.alpha_Y = 0.333;
	p.drs_N = 0;
	p.alpha_N = 0.333;
	p.depreciation = 0.1 / 4.0;
	p.labtax = 0.25;
	p.lumptransfer = 0.05;
	p.corptax = 0;

	// ParamConfig
	p.elast = 6.0;
	p.drs_Y = 0.934;
	p.drs_N = 0.73;
	p.alpha_Y = 0.32;
	p.alpha_N = 0;
	p.depreciation = 0.07 / 4.0;
}

void set_to_fortran_params(SSCalibrator& cal) {
	cal.calibrateLaborDisutility = true;
	cal.calibrateDiscountRate = true;
	cal.calibrateRb = true;
}

int main () {
	std::string income_dir = "2point_3_5";

	Options options; 
	options.fast = false;
	options.print_diagnostics = false;

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

	set_to_fortran_params(params);

	params.setup(options);
	global_params_ptr = &params;

	// Calibration
	SSCalibrator calibrator;
	calibrator.calibrateLaborDisutility = false;
	calibrator.calibrateRb = true;
	calibrator.calibrateDiscountRate = true;
	set_to_fortran_params(calibrator);
	calibrator.setup(params);
	global_calibrator_ptr = &calibrator;

	Model model = Model(params, income_dir);

	// guess rho, chi,labor_occ, capital, and rb
	int n = calibrator.nmoments;
	hank_float_type x[n];
	calibrator.fill_xguess(params, model, x);

	hank_float_type fvec[n];
	double tol = 1.0e-9;

	int lwa = n * (3 * n + 13);
	hank_float_type wa[lwa];

	void *z = NULL;
	cminpack_hybrd1_fnname(fcn, z, n, x, fvec, tol, wa, lwa);

	// int iflag=0;
	// find_initial_steady_state(n, x, fvec, iflag);
}