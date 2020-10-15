#ifndef _INITIAL_STEADY_STATE_H
#define _INITIAL_STEADY_STATE_H

#include <hank_config.h>
#include <vector>
#include <equilibrium.h>

// Forward declarations
class Model;

class Parameters;

// Provides attributes and methods for steady state computations
class SteadyStateBase {
	public:
		enum class SSType { initial, final };

		SSType mode = SSType::initial;

		double output = 1.0;

		double varieties = 1.0;

		double rho, capital, capital_Y, capital_N, price_W;

		double investment, labor_Y, labor_N, labor;

		std::vector<hank_float_type> labor_occ, labshareY, labshareN, labfracY, labfracN;

		hank_float_type capshareY, capshareN, capfracY, capfracN;

		double netprofit_W, grossprofit_R, netprofit_R, profit;

		double rcapital, tfp_N, tfp_Y;

		std::vector<double> wage_occ;

		double wage_Y, wage_N, mc_Y, mc_N;

		double ra, dividend_A, dividend_B, equity_A, equity_B;

		std::vector<double> netwagegrid;

		double taxrev, Enetwage;

		double illprice, illpricedot, illshares;

		double riskaver;

		double mpshock = 0;
};

class SteadyState : public SteadyStateBase {
	public:
		SteadyState(const Parameters& p_, const Model& model_, SSType mode_);

		SteadyState(const SteadyState& other_ss);

		SteadyState& operator=(const SteadyState& other_ss);

		SteadyState& operator=(SteadyState&& other_ss);

		void guess_labor_occ();

		void compute();

		void compute_profits();

		void compute_factor_prices();

		void compute_dividends();

		void compute_govt();

		void print_variables() const;

		const Model& model;

		const Parameters& p;

		Equilibrium equilibrium;
};

#endif