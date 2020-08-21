#ifndef _BELLMAN_H
#define _BELLMAN_H

#include <algorithm>
#include <utilities.h>
#include <hank.h>

struct ValueFnDerivatives {
	static const int StationaryPtOrLimit = -999.9;
	double VaF, VaB, VbF, VbB;
};

namespace Upwinding {
	class ConUpwind;
	class Policies;
}

class Model;
class SteadyState;
class Parameters;

class HJB {
	public:
		HJB(const Model& model_, const SteadyState& ss);

		void iterate(const SteadyState& ss);

		Upwinding::Policies update_policies(const SteadyState& ss);

		void update_value_fn(const SteadyState& ss, const Upwinding::Policies& policies);

		ValueFnDerivatives compute_derivatives(int ia, int ib, int iy) const;

		Upwinding::ConUpwind optimal_consumption(double Vb, double bdrift, double netwage, double chi, double idioscale) const;
		Upwinding::ConUpwind optimal_consumption_no_laborsupply(double Vb, double bdrift, double netwage) const;
		Upwinding::ConUpwind optimal_consumption_sep_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const;
		Upwinding::ConUpwind optimal_consumption_ghh_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const;

		const Model& model;

		const Parameters& p;

		boost3d V;

		int maxiter = 500;
		int dispfreq = 50;
		double vtol = 1.0e-8;
		double delta = 1.0e6;

		double dVamin = 1.0e-8;
		double dVbmin = 1.0e-8;
};

#endif