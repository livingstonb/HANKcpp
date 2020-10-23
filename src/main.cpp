
#include <iostream>
#include <hank_config.h>
#include <parameters.h>
#include <hank_types.h>
#include <model.h>
#include <equilibrium.h>
#include <bellman.h>
#include <stationary_dist.h>
#include <distribution_statistics.h>
#include <adjustment_costs.h>
#include <utilities.h>
#include <math.h>
#include <calibration.h>
#include <impulse_responses.h>
#include <memory>

#include <cminpack_wrapper.h>

const Options *global_hank_options = NULL;

std::string income_dir = "2point_3_5";

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
	p.make_profit_correction = true;

	// ParamConfig
	// p.elast = 6.0;
	// p.drs_Y = 0.934;
	// p.drs_N = 0.73;
	// p.alpha_Y = 0.32;
	// p.alpha_N = 0;
	// p.depreciation = 0.07 / 4.0;

	// ParamConfig 2
	p.elast = 6.89655172413793;
	p.drs_Y = 1.0;
	p.drs_N = 0.310344827586207;
	p.alpha_N = 0.3333333333333;
	p.alpha_Y = 0.3333333333333;
}

void set_to_fortran_params(HANKCalibration::SSCalibrator& cal) {
	cal.calibrateLaborDisutility = true;
	cal.calibrateDiscountRate = true;
	cal.calibrateRb = true;
}

void compute_irfs(const HANKCalibration::ObjectPointers& object_ptrs) {
	const Parameters& p = *object_ptrs.ptr1;
	const Model& model = *object_ptrs.ptr2;
	const EquilibriumInitial& iss = *object_ptrs.ptr3;

	IRF irf(p, model, iss);
	irf.shock.type = ShockType::tfp_Y;
	irf.permanentShock = false;

	irf.setup();

	irf.compute();
}

int main () {
	Options options; 
	options.fast = false;
	options.print_diagnostics = false;
	options.skip_calibration = true;

	global_hank_options = &options;

	HANKCalibration::ObjectPointers object_ptrs;

	// Parameters
	object_ptrs.ptr1.reset(new Parameters);
	Parameters& params = *object_ptrs.ptr1;

	params.income_dir = "2point_3_5";
	params.rho = 0.022;
	params.drs_N = 0;
	params.drs_Y = 0.9;
	params.dmax = 1e6;
	params.borrowing = true;
	// params.deathrate = 0.0;
	params.amax = 400;
	params.na = 25;
	params.nb_pos = 25;
	params.depreciation = 0.05 / 4;
	params.elast = 2;
	params.nocc = 1;
	// params.riskaver = 1.0;
	params.rb = 0.03 / 4.0;

	set_to_fortran_params(params);
	params.setup(options);

	object_ptrs.ptr2.reset(new Model(params));
	const Model& model = *object_ptrs.ptr2;

	if ( options.skip_calibration ) {
		object_ptrs.ptr3.reset(new EquilibriumInitial);
		EquilibriumInitial& iss = *object_ptrs.ptr3;
		iss.solve(params, model);

		HJB hjb(model, iss);
		hjb.iterate(iss);

		StationaryDist sdist;
		sdist.gtol = 1.0e-9;
		sdist.compute(model, iss, hjb);

		object_ptrs.ptr4.reset(new DistributionStatistics(params, model, hjb, sdist));
		const DistributionStatistics& stats = *object_ptrs.ptr4;
		stats.print();

		iss.update_with_stats(stats);

		compute_irfs(object_ptrs);
	}
	else {
		// Calibrate
		object_ptrs.ptr5.reset(new HANKCalibration::SSCalibrator);
		HANKCalibration::SSCalibrator& cal = *object_ptrs.ptr5;
		set_to_fortran_params(cal);
		cal.setup(params);

		// Guess rho, chi,labor_occ, capital, and rb
		int n = cal.nmoments;
		hank_float_type x[n];
		cal.fill_xguess(params, model, x);

		cminpack_hybrd1_wrapper(HANKCalibration::initial_steady_state_obj_fn, &object_ptrs, n, x);
	}

	// compute_irfs(object_ptrs);
}