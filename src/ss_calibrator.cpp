#include <ss_calibrator.h>
#include <functional>
#include <iostream>

using namespace std::placeholders;

namespace {
	double deviation_median_illiq_wealth(const SSCalibrationArgs& args) {
		return args.stats->a_pctiles[5] / args.p->illiqWealthTarget.value - 1.0;
	}

	double deviation_mean_liq_wealth(const SSCalibrationArgs& args) {
		return args.stats->Eb / args.p->liqWealthTarget.value - 1.0;
	}

	double deviation_median_liq_wealth(const SSCalibrationArgs& args) {
		return args.stats->b_pctiles[5] / args.p->liqWealthTarget.value - 1.0;
	}

	double illiq_market_clearing(const SSCalibrationArgs& args) {
		return args.stats->Ea / (args.iss->capital + args.iss->equity_A) - 1.0;
	}

	double labor_market_clearing(const SSCalibrationArgs& args, int io) {
		return args.stats->Elabor_occ[io] / args.iss->labor_occ[io] - 1.0;
	}

	double hours_target(const SSCalibrationArgs& args) {
		return args.stats->Ehours / args.p->hourtarget - 1.0;
	}
}

SSCalibrator::SSCalibrator(const Parameters &p) {
	if ( p.illiqWealthTarget.is_median() )
		obj_functions.push_back(deviation_median_illiq_wealth);

	if ( p.liqWealthTarget.is_mean() )
		obj_functions.push_back(deviation_mean_liq_wealth);
	else if ( p.liqWealthTarget.is_median() )
		obj_functions.push_back(deviation_median_liq_wealth);

	// Market clearing for illiquid asset
	obj_functions.push_back(illiq_market_clearing);

	// Labor market clearing
	for (int io=0; io<p.nocc; ++io) {
		deviation_fn_type labclearing_io = std::bind(labor_market_clearing, _1, io);
		obj_functions.push_back(labclearing_io);
	}

	// Hours target
	obj_functions.push_back(hours_target);
}

void SSCalibrator::fill_fvec(const SSCalibrationArgs& args, double fvec[]) const {
	for (int i=0; i<nmoments(); ++i) {
		fvec[i] = obj_functions[i](args);
	}
}

void SSCalibrator::fill_xguess(const Parameters &p, const Model& model, double xvec[]) {
	int ix = 0;

	// Labor market clearing
	for (int io=0; io<p.nocc; ++io) {
		xvec[ix] = log(p.hourtarget * p.meanlabeff * model.occdist[io]);
		ix_occdist.push_back(ix);
		++ix;
	}

	// Discount rate
	if ( p.calibrateDiscountRate ) {
		check_size(ix);
		xvec[ix] = log(p.rho);
		ix_rho = ix;
		++ix;
	}

	// Liquid returns
	if ( p.calibrateRb ) {
		check_size(ix);
		xvec[ix] = log(p.rb);
		ix_rb = ix;
		++ix;
	}

	// Labor disutility
	if ( p.calibrateLaborDisutility ) {
		check_size(ix);
		xvec[ix] = log(p.chi);
		ix_chi = ix;
		++ix;
	}

	if ( ix == nmoments() - 1) {
		// Use capital as the last moment
		xvec[ix] = log(p.target_KY_ratio);
		ix_capital = ix;
		++ix;
	}
	else if ( ix < nmoments() - 1 ) {
		std::cerr << "Too few moments\n";
		throw 0;
	}
	else if ( ix > nmoments() ) {
		std::cerr << "Too many moments\n";
		throw 0;
	}
}

void SSCalibrator::check_size(int ix) const {
	if ( ix >= nmoments() ) {
		std::cerr << "Too many guesses\n";
		throw 0;
	}
}

void SSCalibrator::update_params(Parameters *p, double xvec[]) const {
	if ( ix_rho > 0 ) {
		p->rho = exp(xvec[ix_rho]);
		std::cout << "  rho = " << p->rho << '\n';
	}

	if ( ix_rb > 0 ) {
		p->rb = exp(xvec[ix_rb]);
		std::cout << "  rb = " << p->rb << '\n';
	}

	if ( ix_chi > 0 ) {
		p->chi = exp(xvec[ix_chi]);
		std::cout << "  chi = " << p->chi << '\n';
	}

	p->update();
}

void SSCalibrator::update_ss(const Parameters& p, SteadyState *iss, double xvec[]) const {
	for (int io=0; io<ix_occdist.size(); ++io) {
		iss->labor_occ.push_back(exp(xvec[ix_occdist[io]]));
		std::cout << "  labor_" << io << " = " << exp(xvec[ix_occdist[io]]) << '\n';
	}

	if ( ix_capital > 0 ) {
		iss->capital = exp(xvec[ix_capital]);
		std::cout << "  capital = " << iss->capital << '\n';
	}
	else
		iss->capital = p.target_KY_ratio;
}