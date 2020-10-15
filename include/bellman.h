#ifndef _BELLMAN_H
#define _BELLMAN_H

#include <hank_config.h>
#include <hank_types.h>
#include <upwinding.h>
#include <math.h>

// Forward declarations
class Model;

class EquilibriumElement;

class Parameters;

// Container for value function derivatives
struct ValueFnDerivatives {
	static const int StationaryPtOrLimit = -999.9;

	double VaF, VaB, VbF, VbB;
};

namespace Bellman {
	class Drifts {
		public:
			Drifts() {}

			Drifts(double s, double d, double areturn, double acost, bool kfe);

			double aB, aF, bB, bF;
	};
}

// Class for solving the HJB
class HJB {
	private:
		Upwinding::Policies update_policies(const EquilibriumElement& ss);

		void update_value_fn(const EquilibriumElement& ss, const Upwinding::Policies& policies);

		ValueFnDerivatives compute_derivatives(int ia, int ib, int iy) const;

		Upwinding::ConUpwind optimal_consumption(double Vb, double bdrift, double netwage, double chi, double idioscale) const;

		Upwinding::ConUpwind optimal_consumption_no_laborsupply(double Vb, double bdrift, double netwage) const;

		Upwinding::ConUpwind optimal_consumption_sep_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const;

	public:
		HJB(const Model& model_, const EquilibriumElement& ss);

		void iterate(const EquilibriumElement& ss);

		void print_variables() const;

		const Model& model;

		const Parameters& p;

		vector3dr V;

		int maxiter = 500;

		int dispfreq = 50;

		double vtol = 1.0e-8;

		double delta = 1.0e6;

		double dVamin = 1.0e-8;

		double dVbmin = 1.0e-8;

		double riskaver;

		Upwinding::Policies optimal_decisions;
};

#endif