#ifndef _INITIAL_STEADY_STATE_H
#define _INITIAL_STEADY_STATE_H

#include <iostream>
#include <mkl.h>
#include <model.h>
#include <parameters.h>

void solve_initial_steady_state(const Model& model);

void solve_initial_prices(const Model& model);

double compute_ss_capital_output_ratio(const Parameters& params, double lprice_W);

#endif