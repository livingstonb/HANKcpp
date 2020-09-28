#ifndef _IMPULSE_RESPONSES_H
#define _IMPULSE_RESPONSES_H

#include <hank_types.h>
#include <hank_eigen_dense.h>

class Parameters;

class Model;

class SteadyState;

class TransEquilibrium {
	public:
		VectorXr tfp_Y, mpshock;
};

enum class ShockType { tfp_Y, monetary, riskaver, none };

class TransShock {
	private:
		const double value_not_set = -999.9;

	public:
		ShockType type = ShockType::none;

		double size = value_not_set;

		double pers = value_not_set;

		void setup();
};

class IRF {
	public:
		IRF(const Parameters& p_, const Model& model_, const SteadyState& iss_);

		void setup();

		void construct_delta_trans_vectors();

		void compute();

		TransShock shock;

		double shock_size = 0;

		int Ttrans = 30;

		int nendtrans = 5; // Number of steps to enforce return to steady state

		int npricetrans;

		double deltatransmin = 1.0 / 3.0;

		double deltatranstot = 400;

		bool solveFlexPriceTransitions = false;

		bool solveStickyPriceTransitions = true;

		VectorXr deltatransvec, cumdeltatrans;

		TransEquilibrium trans_equm;

		const Parameters& p;

		const Model& model;

		const SteadyState& iss;
};


#endif