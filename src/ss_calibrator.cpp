#include <ss_calibrator.h>
#include <functional>

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
	for (int i=0; i<nvals(); ++i) {
		fvec[i] = obj_functions[i](args);
	}
}