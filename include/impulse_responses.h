#ifndef _IMPULSE_RESPONSES_H
#define _IMPULSE_RESPONSES_H

#include <hank_types.h>
#include <hank_eigen_dense.h>
#include <memory>
#include <equilibrium.h>

class Parameters;

class Model;

enum class ShockType { tfp_Y, monetary, riskaver, none };

class TransShock {
	public:
		ShockType type = ShockType::none;

		hank_float_type size = HANK::ValueNotSet;

		hank_float_type pers = HANK::ValueNotSet;

		void setup();
};

class IRF {
	public:
		IRF(const Parameters& p_, const Model& model_, const EquilibriumElement& iss_);

		enum class SolverType { hybrd1, broyden };

		void setup();

		void construct_delta_trans_vectors();

		void compute();

		void transition_fcn(int n, const hank_float_type *x, hank_float_type *z);

		void make_transition_guesses(int n, const hank_float_type *x, hank_float_type *z);

		void set_shock_paths();

		SolverType solver = SolverType::broyden;

		TransShock shock;

		double shock_size = 0;

		int Ttrans = 30;

		int nendtrans = 5; // Number of steps to enforce return to steady state

		int npricetrans;

		double deltatransmin = 1.0 / 3.0;

		double deltatranstot = 400;

		bool solveFlexPriceTransitions = false;

		bool solveStickyPriceTransitions = true;

		bool permanentShock = false;

		double min_price_W = 0.05;

		double max_price_W = 1.0 - 1.0e-6;

		VectorXr deltatransvec, cumdeltatrans;

		TransEquilibrium trans_equm;

		EquilibriumElement final_equm;

		const Parameters& p;

		const Model& model;

		const EquilibriumElement& initial_equm;
};


#endif