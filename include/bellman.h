#ifndef _BELLMAN_H
#define _BELLMAN_H

#include <model.h>
#include <steady_state.h>

class HJB {
	public:
		HJB(const Model& model, const SteadyState& ss);

		void iterate(const Model& model, const SteadyState& ss);

		boost_array_type<double, 3> V;
};

boost_array_type<double, 3> make_value_guess(const Model& model, const SteadyState& ss);



#endif