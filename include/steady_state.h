#ifndef _STEADY_STATE_H
#define _STEADY_STATE_H

#include <iostream>
#include <mkl.h>
#include <model.h>
#include <parameters.h>
#include <options.h>

void solve_initial_steady_state(const Model& model, const Options& options);

void solve_initial_prices(const Model& model);

double compute_ss_capital_output_ratio(const Model& model, double lprice_W);

class Prices {
	public:
		double ra;
};

#endif