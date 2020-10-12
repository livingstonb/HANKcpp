#ifndef _INITIAL_STEADY_STATE_H
#define _INITIAL_STEADY_STATE_H

#include <hank_config.h>
#include <vector>

// Forward declarations
class Model;

class Parameters;

// Provides attributes and methods for steady state computations
class SteadyState {
	public:
		enum class SSType { initial, final };

		SteadyState(const Parameters& p_, const Model& model_, SSType mode_);

		SteadyState& operator=(const SteadyState& other_ss) {
			*this = other_ss;
			return *this;
		}

		void guess_labor_occ();

		void compute();

		void compute_profits();

		void compute_factor_prices();

		void compute_dividends();

		void compute_govt();

		void print_variables() const;

		const Model& model;

		const Parameters& p;

		SSType mode = SSType::initial;

		double output = 1.0;

		double varieties = 1.0;

		double rho, capital, capital_Y, capital_N, price_W;

		double investment, labor_Y, labor_N, labor;

		std::vector<double> labor_occ, labshareY, labshareN, labfracY, labfracN;

		double capshareY, capshareN, capfracY, capfracN;

		double netprofit_W, grossprofit_R, netprofit_R, profit;

		double rcapital, tfp_N, tfp_Y;

		std::vector<double> wage_occ;

		double wage_Y, wage_N, mc_Y, mc_N;

		double ra, dividend_A, dividend_B, equity_A, equity_B;

		std::vector<double> netwagegrid;

		double taxrev, Enetwage;

		double illprice, illpricedot, illshares;

		double chi, riskaver;

		double mpshock = 0;
};

#endif