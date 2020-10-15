#ifndef _IMPULSE_RESPONSES_H
#define _IMPULSE_RESPONSES_H

#include <hank_types.h>
#include <hank_eigen_dense.h>

class Parameters;

class Model;

class SteadyState;

class TransEquilibrium {
	public:
		TransEquilibrium() {}

		void set_array_sizes(const Parameters& p, int T);

		VectorXr tfp_Y, mpshock, riskaver, rb, pi, qcapital, rnom, pricelev;

		VectorXr output;

		MatrixXr labor_occ;

		double ss_riskaver;
};

class Equilibrium {
	public:
		Equilibrium() {}

		Equilibrium(const Parameters& p, const SteadyState& ss);

		double rb, pi, output, qcapital, rnom;

		std::vector<double> labor_occ;
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

		enum class SolverType { hybrd1, broyden };

		void setup();

		void construct_delta_trans_vectors();

		void compute();

		void transition_fcn(int n, const hank_float_type *x, hank_float_type *z);

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

		double trans_riskaver;

		VectorXr deltatransvec, cumdeltatrans;

		TransEquilibrium trans_equm;

		Equilibrium initial_equm, final_equm;

		const Parameters& p;

		const Model& model;

		const SteadyState& iss;
};


#endif