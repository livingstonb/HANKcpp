#ifndef _IMPULSE_RESPONSES_H
#define _IMPULSE_RESPONSES_H

#include <hank_config.h>
#include <hank.h>
#include <memory>
#include <equilibrium.h>

class Parameters;

class Model;

class IRF;

enum class ShockType { tfp_Y, monetary, riskaver, none };

class TransShock {
	private:
		void setup();

	public:
		ShockType type = ShockType::none;

		hank_float_type size = HANK::ValueNotSet;

		hank_float_type pers = HANK::ValueNotSet;

		friend class IRF;
};

class IRF {
	private:
		const Parameters& p;

		const Model& model;

		const EquilibriumInitial& initial_equm;

		void transition_fcn(int n, const hank_float_type *x, hank_float_type *z);

		void compute_remaining_variables();

		void set_shock_paths();

		void find_final_steady_state();

	public:
		IRF(const Parameters& p_, const Model& model_, const EquilibriumInitial& iss_);

		enum class SolverType { hybrd1, broyden };

		void setup();

		void compute();

		SolverType solver = SolverType::broyden;

		TransShock shock;

		double shock_size = 0;

		int Ttrans = 30;

		int nendtrans = 5; // Number of steps to enforce return to steady state

		int npricetrans;

		double deltatransmin = 1.0 / 3.0;

		double deltatranstot = 400;

		bool flextransition = false;

		bool stickytransition = true;

		bool permanentShock = false;

		double min_price_W = 0.05;

		double max_price_W = 1.0 - 1.0e-6;

		std::vector<hank_float_type> deltatransvec, cumdeltatrans;

		std::vector<EquilibriumTrans> trans_equm;

		std::shared_ptr<EquilibriumFinal> final_equm_ptr = nullptr;
};

#endif