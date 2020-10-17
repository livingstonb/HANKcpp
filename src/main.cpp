
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
#include <ss_calibrator.h>
#include <impulse_responses.h>
#include <memory>

#include <cminpack_wrapper.h>

const Options *global_hank_options = NULL;

std::string income_dir = "2point_3_5";

std::unique_ptr<SSCalibrator> global_calibrator_ptr(nullptr);

using HANKObjectPointers = UniquePtrContainer<Parameters, Model, EquilibriumElement, DistributionStatistics, hank_float_type[]>;

int fcn(void* args_void_ptr, int n, const real *x, real *fvec, int /* iflag */ ) {
	assert(global_calibrator_ptr->nmoments == n);

	std::cout << "\nCalibration parameters updated:\n";

	HANKObjectPointers& args = *(HANKObjectPointers *) args_void_ptr;
	Parameters& p = *args.ptr1;

	global_calibrator_ptr->update_params(&p, x);

	args.ptr2 = std::make_unique<Model>(p, income_dir);
	Model& model = *args.ptr2;

	args.ptr3 = std::make_unique<EquilibriumElement>();
	EquilibriumElement& iss = *args.ptr3;
	global_calibrator_ptr->update_ss(&p, &iss, x);
	std::cout << '\n';

	iss.create_initial_steady_state(p, model);

	HJB hjb(*args.ptr2, iss);
	hjb.iterate(iss);

	StationaryDist sdist;
	sdist.gtol = 1.0e-9;
	sdist.compute(model, iss, hjb);

	args.ptr4 = std::make_unique<DistributionStatistics>(p, model, hjb, sdist);
	DistributionStatistics& stats = *args.ptr4;
	stats.print();

	SSCalibrationArgs cal_args(&p, &model, &stats, &iss);
	global_calibrator_ptr->fill_fvec(cal_args, fvec);
	global_calibrator_ptr->print_fvec(fvec);

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
	options.skip_calibration = true;

	global_hank_options = &options;

	HANKObjectPointers object_ptrs;

	// Parameters
	object_ptrs.ptr1 = std::make_unique<Parameters>();
	Parameters& params = *object_ptrs.ptr1;

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

	object_ptrs.ptr2 = std::make_unique<Model>(params, income_dir);
	Model& model = *object_ptrs.ptr2;
	// global_current_model_ptr.reset(new Model(params, income_dir));
	// Model& model = *global_current_model_ptr;

	// global_current_iss_ptr.reset(new Equilibrium(1));

	if ( options.skip_calibration ) {
		object_ptrs.ptr3 = std::make_unique<EquilibriumElement>();
		EquilibriumElement& iss = *object_ptrs.ptr3;
		iss.create_initial_steady_state(params, model);

		HJB hjb(model, iss);
		hjb.iterate(iss);

		StationaryDist sdist;
		sdist.gtol = 1.0e-9;
		sdist.compute(model, iss, hjb);

		object_ptrs.ptr4 = std::make_unique<DistributionStatistics>(params, model, hjb, sdist);
		object_ptrs.ptr4->print();
	}
	else {
		// Calibrate
		global_calibrator_ptr.reset(new SSCalibrator);
		global_calibrator_ptr->calibrateLaborDisutility = false;
		global_calibrator_ptr->calibrateRb = true;
		global_calibrator_ptr->calibrateDiscountRate = true;
		set_to_fortran_params(*global_calibrator_ptr);
		global_calibrator_ptr->setup(params);

		// Guess rho, chi,labor_occ, capital, and rb
		int n = global_calibrator_ptr->nmoments;
		hank_float_type x[n];
		global_calibrator_ptr->fill_xguess(params, model, x);

		int info = cminpack_hybrd1_wrapper(fcn, (void *) &object_ptrs, n, x);
		HankUtilities::check_cminpack_success(info);
	}

	EquilibriumElement& iss = *object_ptrs.ptr3;
	IRF irf(params, *object_ptrs.ptr2, iss);
	irf.shock.type = ShockType::tfp_Y;
	irf.permanentShock = true;
	irf.setup();

	irf.compute();
}