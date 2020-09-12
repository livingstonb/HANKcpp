#ifndef _INITIAL_STEADY_STATE_H
#define _INITIAL_STEADY_STATE_H

#include <hank_config.h>
#include <vector>

// Forward declarations
class Model;

class Parameters;

// Provides attributes and methods for steady state computations
class SteadyState {
	private:
		const Model& model;
		const Parameters& p;

	public:
		SteadyState(const Model& model_);

		void update();

		// Normalize to 1
		double output = 1.0;
		double varieties = 1.0;
		double totoutput;

		// Profits, prices, and dividends
		double price_W, grossprofit_W, netprofit_W, grossprofit_R, netprofit_R;
		double profit, dividend_A, dividend_B, equity_A, equity_B;

		// Production
		double tfp_N, tfp_Y, mc_N, mc_Y;

		// Labor market
		double labor_Y, labor_N, wage_N, wage_Y;
		std::vector<double> netwagegrid;

		// Capital
		double ra, rcapital, capital, capital_N, capital_Y;
		double K_totoutput_ratio;

		// Labor disutility
		double chi;

		// Productivity grid
		std::vector<double> yprodgrid;
};




#endif