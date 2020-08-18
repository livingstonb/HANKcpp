#ifndef _BELLMAN_H
#define _BELLMAN_H

#include <algorithm>
#include <HankNumerics.h>

#include <model.h>
#include <steady_state.h>

void check_progress(double vdiff, int freq, int i, double vtol);

boost_array_type<double, 3> make_value_guess(const Model& model, const SteadyState& ss);

struct ConUpwind {
	double c, h, s, Hc;
};

const int STATIONARY_PT_OR_LIMIT = -999.9;

constexpr bool is_stationary_pt_or_limit(double Vb) {
	return (Vb <= STATIONARY_PT_OR_LIMIT);
}


class HJB {
	public:
		HJB(const Model& model_, const SteadyState& ss);

		void iterate(const SteadyState& ss);

		void update(const SteadyState& ss);

		void compute_derivatives(
			double& VaF, double& VbF, double& VaB, double& VbB,
			int ia, int ib, int iy) const;

		ConUpwind optimal_consumption(double Vb, double bdrift, double netwage, double chi, double idioscale) const;
		ConUpwind optimal_consumption_no_laborsupply(double Vb, double bdrift, double netwage) const;
		ConUpwind optimal_consumption_sep_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const;
		ConUpwind optimal_consumption_ghh_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const;

		const Model& model;

		const Parameters& p;

		boost_array_type<double, 3> V;

		int maxiter = 500;
		int dispfreq = 50;
		double vtol = 1.0e-8;

		double dVamin = 1.0e-8;
		double dVbmin = 1.0e-8;
};




#endif