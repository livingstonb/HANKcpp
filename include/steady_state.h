#ifndef _INITIAL_STEADY_STATE_H
#define _INITIAL_STEADY_STATE_H

#include <iostream>
#include <mkl.h>
#include <model.h>
#include <parameters.h>
#include <options.h>

class InitialSteadyState {
	public:
		double ra;
		double profit;

		double dividend_A;
		double dividend_B;
		double equity_A;
		double equity_B;

		double_vector netwagegrid;

		double wage_N;
		double wage_Y;

		double chi;
};

InitialSteadyState find_steady_state(const Model& model, const Options& options);

/*
	Solves for various steady state quantities and returns and object with the
	results.
*/
InitialSteadyState initialize_steady_state(const Model& model);

/*
	Computes the steady state capital-output ratio K/NY, given the model
	parameters and the target for mean illiquid wealth.
*/
double compute_ss_capital_output_ratio(double price_W,
	double targetMeanIll, double depreciation, double alpha_Y,
	double drs_Y, double alpha_N, double drs_N);


#endif