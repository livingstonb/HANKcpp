#ifndef _INITIAL_STEADY_STATE_H
#define _INITIAL_STEADY_STATE_H

#include <iostream>
// #include <mkl.h>
#include <model.h>
#include <parameters.h>
#include <options.h>

double marginal_cost(double tfp, double r, double alpha, double wage);

double compute_ss_capital_output_ratio(double price_W,
	double targetMeanIll, double depreciation, double alpha_Y,
	double drs_Y, double alpha_N, double drs_N);

double compute_ss_capital_output_ratio(const Parameters& p, double price_w);

double quadratic_formula(double a, double b, double c);

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
		double price_W;
		double grossprofit_W, netprofit_W;
		double grossprofit_R, netprofit_R;
		double profit;
		double dividend_A, dividend_B, equity_A, equity_B;

		// Production
		double tfp_N, tfp_Y, mc_N, mc_Y;

		// Labor market
		double labor_Y, labor_N, wage_N, wage_Y;
		double_vector netwagegrid;

		// Capital
		double ra, rcapital;
		double capital, capital_N, capital_Y;
		double K_totoutput_ratio;

		// Labor disutility
		double chi;

		// Productivity grid
		double_vector yprodgrid;
};




#endif