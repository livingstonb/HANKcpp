#ifndef _BELLMAN_H
#define _BELLMAN_H

#include <algorithm>

#include <model.h>
#include <steady_state.h>

void check_progress(double vdiff, int freq, int i, double vtol);

boost_array_type<double, 3> make_value_guess(const Model& model, const SteadyState& ss);

class HJB {
	public:
		HJB(const Model& model_, const SteadyState& ss);

		void iterate(const SteadyState& ss);

		void update(const SteadyState& ss);

		const Model& model;

		boost_array_type<double, 3> V;

		int maxiter = 500;
		int dispfreq = 50;
		double vtol = 1.0e-8;

		double dVamin = 1.0e-8;
		double dVbmin = 1.0e-8;
};




#endif