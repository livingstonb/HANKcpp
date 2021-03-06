#include <calibration.h>
#include <iostream>
#include <string>
#include <math.h>
#include <hank_config.h>
#include <hank.h>

#include <parameters.h>
#include <model.h>
#include <bellman.h>
#include <distribution_statistics.h>
#include <stationary_dist.h>
#include <upwinding.h>
#include <hank_numerics_cminpack.h>

#include <assert.h>

using namespace std::placeholders;

using namespace HANKCalibration;

namespace {
	void check_size(int nmoments, int ix);

	void perform_calibrator_assertions(const SSCalibrator* cal);

	double deviation_median_illiq_wealth(const CalibrationArgs& args);

	double value_median_illiq_wealth(const CalibrationArgs& args);

	double deviation_mean_liq_wealth(const CalibrationArgs& args);

	double value_mean_liq_wealth(const CalibrationArgs& args);

	double deviation_median_liq_wealth(const CalibrationArgs& args);

	double value_median_liq_wealth(const CalibrationArgs& args);

	double illiq_market_clearing(const CalibrationArgs& args);

	double value_mean_illiq_wealth(const CalibrationArgs& args);

	double labor_market_clearing(const CalibrationArgs& args, int io);

	double value_labor(const CalibrationArgs& args, int io);

	double hours_target(const CalibrationArgs& args);

	double value_hours(const CalibrationArgs& args);
}

void SSCalibrator::setup(const Parameters &p)
{
	// Market clearing conditions
	obj_functions.push_back(illiq_market_clearing);
	moment_descriptions.push_back("Capital market clearing");
	variable_names.push_back("E[a]");
	variable_values.push_back(value_mean_illiq_wealth);

	for (int io=0; io<p.nocc; ++io) {
		deviation_fn_type labclearing_io = std::bind(labor_market_clearing, _1, io);
		obj_functions.push_back(labclearing_io);
		moment_descriptions.push_back("Labor market clearing, occ_" + std::to_string(io) + "");
		variable_names.push_back("E[labor_" + std::to_string(io) + "]");

		deviation_fn_type labclearing_val_io = std::bind(value_labor, _1, io);
		variable_values.push_back(labclearing_val_io);
	}

	if ( p.illiqWealthTarget.is_median() ) {
		obj_functions.push_back(deviation_median_illiq_wealth);
		moment_descriptions.push_back("Median illiq wealth");
		variable_names.push_back("Median(a)");
		variable_values.push_back(value_median_illiq_wealth);
	}

	if ( p.liqWealthTarget.is_mean() ) {
		obj_functions.push_back(deviation_mean_liq_wealth);
		moment_descriptions.push_back("Mean liq wealth");
		variable_names.push_back("E[b]");
		variable_values.push_back(value_mean_liq_wealth);
	}
	else if ( p.liqWealthTarget.is_median() ) {
		obj_functions.push_back(deviation_median_liq_wealth);
		moment_descriptions.push_back("Median liq wealth");
		variable_names.push_back("Median(b)");
		variable_values.push_back(value_median_liq_wealth);
	}

	// Hours target
	if ( calibrateLaborDisutility ) {
		obj_functions.push_back(hours_target);
		moment_descriptions.push_back("Hours worked");
		variable_names.push_back("E[hours");
		variable_values.push_back(value_hours);
	}

	// Set number of moments
	nmoments = obj_functions.size();

	// perform_calibrator_assertions(this);
}

void SSCalibrator::fill_fvec(const CalibrationArgs& args, hank_float_type fvec[]) const
{
	for (int i=0; i<nmoments; ++i) {
		fvec[i] = obj_functions[i](args);
	}
}

void SSCalibrator::fill_xguess(const Parameters &p, const Model& model, hank_float_type xvec[])
{
	int ix = 0;

	// Labor inputs
	for (int io=0; io<p.nocc; ++io) {
		xvec[ix] = p.hourtarget * p.meanlabeff * model.occdist[io];
		ix_labor_occ.push_back(ix);
		++ix;
	}

	// Discount rate
	if ( calibrateDiscountRate ) {
		check_size(nmoments, ix);
		xvec[ix] = log(p.rho);
		ix_rho = ix;
		++ix;
	}

	// Liquid returns
	if ( calibrateRb ) {
		check_size(nmoments, ix);
		xvec[ix] = log(p.rb);
		ix_rb = ix;
		++ix;
	}

	// Labor disutility
	if ( calibrateLaborDisutility ) {
		check_size(nmoments, ix);
		xvec[ix] = p.chi;
		ix_chi = ix;
		++ix;
	}

	if ( ix == nmoments - 1) {
		// Use capital as the last input variable
		xvec[ix] = p.target_KY_ratio;
		ix_capital = ix;
	}
	else if ( ix < nmoments - 1 ) {
		std::cerr << "Too few moments\n";
		throw 0;
	}
	else if ( ix > nmoments ) {
		std::cerr << "Too many moments\n";
		throw 0;
	}
}

void SSCalibrator::update_params(Parameters *p, const hank_float_type *xvec) const
{
	if ( ix_rho >= 0 ) {
		p->rho = exp(xvec[ix_rho]);
		std::cout << "  rho = " << p->rho << '\n';
	}

	if ( ix_rb >= 0 ) {
		p->rb = exp(xvec[ix_rb]);
		std::cout << "  rb = " << p->rb << '\n';
	}

	if ( ix_chi >= 0 ) {
		p->chi = xvec[ix_chi];
		std::cout << "  chi = " << p->chi << '\n';
	}

	p->update();
}

void SSCalibrator::update_ss(const Parameters* p, EquilibriumInitial *iss, const hank_float_type *xvec) const
{
	for (unsigned int io=0; io<ix_labor_occ.size(); ++io) {
		iss->labor_occ.push_back(xvec[ix_labor_occ[io]]);
		std::cout << "  labor_" << io << " = " << xvec[ix_labor_occ[io]] << '\n';
	}

	if ( ix_capital >= 0 ) {
		iss->capital = xvec[ix_capital];
		std::cout << "  capital = " << iss->capital << '\n';
	}
	else
		iss->capital = p->target_KY_ratio;
}

void SSCalibrator::print_fvec(hank_float_type* fvec) const {
	hank_float_type norm = 0;
	for (unsigned int im=0; im<moment_descriptions.size(); ++im) {
		std::cout << moment_descriptions[im];
		std::cout << "  " << fvec[im] << '\n';
		norm += pow(fvec[im], 2.0);
	}
	norm = sqrt(norm);

	std::cout << "\n fnorm = " << norm << "\n";
	std::cout << "--------------------------\n\n";
}

int HANKCalibration::initial_steady_state_obj_fn(void* args_void_ptr, int n, const hank_float_type *x, hank_float_type *fvec, int /* iflag */ )
{
	std::cout << "\nCalibration parameters updated:\n";

	ObjectPointers& args = *(ObjectPointers *) args_void_ptr;
	Parameters& p = *args.ptr1;
	SSCalibrator& cal = *args.ptr5;

	assert(cal.nmoments == n);

	cal.update_params(&p, x);

	args.ptr2.reset(new Model(p));
	const Model& model = *args.ptr2;

	args.ptr3.reset(new EquilibriumInitial);
	EquilibriumInitial& iss = *args.ptr3;
	cal.update_ss(&p, &iss, x);
	std::cout << '\n';

	iss.setup(p, model);
	iss.solve(p);

	HJB hjb(p, model, iss);
	hjb.iterate(iss);

	StationaryDist sdist;
	sdist.gtol = 1.0e-9;
	sdist.compute(p, model, iss, *hjb.optimal_decisions);

	args.ptr4.reset(new DistributionStatistics(p, model, *hjb.optimal_decisions, sdist.density));
	const DistributionStatistics& stats = *args.ptr4;
	HANK::print(&stats);

	iss.update_after_solving(stats, hjb);

	CalibrationArgs cal_args(args);
	cal.fill_fvec(cal_args, fvec);
	// cal.print_fvec(fvec);

	std::vector<hank_float_type> variable_values;
	for (int i=0; i<n; ++i)
		variable_values.push_back(cal.variable_values[i](cal_args));

	if ( cal.printDetailed ) {
		HANK::OptimStatus optim_status(cal.moment_descriptions, cal.variable_names, fvec, variable_values, cal.iter);
		HANK::print(&optim_status);
	}
	else {
		HANK::OptimNorm optim_norm(fvec, n, cal.iter);
		HANK::print(&optim_norm);
	}

	++cal.iter;

	return 0;
}

namespace {
	void check_size(int nmoments, int ix)
	{
		if ( ix >= nmoments ) {
			std::cerr << "Too many guesses\n";
			throw 0;
		}
	}

	void perform_calibrator_assertions(const SSCalibrator* cal)
	{
		if ( (!cal->calibrateLaborDisutility) & (cal->ix_labor_occ.size() > 0) ) {
			std::cerr << "Labor disutility calibration off but ix_labor_occ is non-empty";
			throw 0;
		}
		else if ( cal->calibrateLaborDisutility & (cal->ix_labor_occ.size() == 0) ) {
			std::cerr << "Labor disutility calibration on but ix_labor_occ is empty";
			throw 0;
		}

		if ( (!cal->calibrateRb) & (cal->ix_rb >= 0) ) {
			std::cerr << "Rb calibration off but ix_rb is non-negative";
			throw 0;
		}
		else if ( cal->calibrateRb & (cal->ix_rb < 0) ) {
			std::cerr << "Rb calibration on but ix_rb is negative";
			throw 0;
		}

		if ( (!cal->calibrateDiscountRate) & (cal->ix_rho >= 0) ) {
			std::cerr << "Discount rate calibration off but ix_rho is non-negative";
			throw 0;
		}
		else if ( cal->calibrateDiscountRate & (cal->ix_rho < 0) ) {
			std::cerr << "Discount rate calibration on but ix_rho is negative";
			throw 0;
		}
	}

	double deviation_median_illiq_wealth(const CalibrationArgs& args)
	{
		return args.stats.a_pctiles[5] / args.p.illiqWealthTarget.value - 1.0;
	}

	double value_median_illiq_wealth(const CalibrationArgs& args)
	{
		return args.stats.a_pctiles[5];
	}

	double deviation_mean_liq_wealth(const CalibrationArgs& args)
	{
		return args.stats.Eb / args.p.liqWealthTarget.value - 1.0;
	}

	double value_mean_liq_wealth(const CalibrationArgs& args)
	{
		return args.stats.Eb;
	}

	double deviation_median_liq_wealth(const CalibrationArgs& args)
	{
		return args.stats.b_pctiles[5] / args.p.liqWealthTarget.value - 1.0;
	}

	double value_median_liq_wealth(const CalibrationArgs& args)
	{
		return args.stats.b_pctiles[5];
	}

	double illiq_market_clearing(const CalibrationArgs& args)
	{
		return args.stats.Ea / (args.iss.capital + args.iss.equity_A) - 1.0;
	}

	double value_mean_illiq_wealth(const CalibrationArgs& args)
	{
		return args.stats.Ea;
	}

	double labor_market_clearing(const CalibrationArgs& args, int io)
	{
		return args.stats.Elabor_occ[io] * args.model.occdist[io] / args.iss.labor_occ[io] - 1.0;
	}

	double value_labor(const CalibrationArgs& args, int io)
	{
		return args.stats.Elabor_occ[io];
	}

	double hours_target(const CalibrationArgs& args)
	{
		return (args.stats.Ehours / args.p.hourtarget - 1.0) / 100.0;
	}

	double value_hours(const CalibrationArgs& args)
	{
		return args.stats.Ehours;
	}
}