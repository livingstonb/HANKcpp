#ifndef _BELLMAN_H
#define _BELLMAN_H

#include <algorithm>
#include <HankNumerics.h>

#include <model.h>
#include <steady_state.h>

struct ConUpwind {
	double c, h, s, Hc;
	bool valid = false;
};

struct DepositUpwind {
	double d, Hd;
	bool valid = false;
};

struct ValueFnDerivatives {
	static const int StationaryPtOrLimit = -999.9;
	double VaF, VaB, VbF, VbB;
};

class HJB {
	public:
		HJB(const Model& model_, const SteadyState& ss);

		struct NoLaborSupply {};
		struct SepLabor {};
		struct GHHLabor {};

		void iterate(const SteadyState& ss);

		void update(const SteadyState& ss);

		ValueFnDerivatives compute_derivatives(int ia, int ib, int iy) const;

		ConUpwind optimal_consumption(double Vb, double bdrift, double netwage, double chi, double idioscale) const;
		ConUpwind optimal_consumption_no_laborsupply(double Vb, double bdrift, double netwage) const;
		ConUpwind optimal_consumption_sep_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const;
		ConUpwind optimal_consumption_ghh_labor(double Vb, double bdrift, double netwage, double chi, double idioscale) const;

		const Model& model;

		const Parameters& p;

		boost3d V;

		int maxiter = 500;
		int dispfreq = 50;
		double vtol = 1.0e-8;

		double dVamin = 1.0e-8;
		double dVbmin = 1.0e-8;
};

#endif