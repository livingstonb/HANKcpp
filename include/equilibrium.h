#ifndef _EQUILIBRIUM_H
#define _EQUILIBRIUM_H

#include <hank_config.h>
#include <hank_types.h>
#include <memory>

class Parameters;

class Model;

class Equilibrium {
	public:
		Equilibrium() {}

		void create_transition();

		virtual void set_from_parameters(const Parameters& p, const Model& model);

		virtual void compute_factors(const Model& model);

		void compute_profits();

		void compute_factor_prices();

		void compute_dividends(const Parameters& p);

		void compute_netwage(const Parameters& p, const Model& model);

		hank_float_type alpha_Y, alpha_N, price_W, drs_Y, drs_N, riskaver, rho;

		hank_float_type capshareY, capshareN, capfracY, capfracN, capital_Y, capital_N;

		hank_float_type valcapital, labor_Y, labor_N, labor;

		hank_float_type tfp_Y, tfp_N, investment, illprice, illshares, illpricedot;

		hank_float_type netprofit_W, grossprofit_R, netprofit_R, profit;

		hank_float_type rcapital, wage_Y, wage_N, mc_Y, mc_N;

		hank_float_type ra, rb, rnom, pi, dividend_A, dividend_B, equity_A, equity_B;

		hank_float_type Enetwage, taxrev, transfershock, lumptransfer, rborr;

		hank_float_type bond, govbond, labtax, govexp;

		hank_float_type output = 1.0;

		hank_float_type varieties = 1.0;

		hank_float_type qcapital = 1.0;

		hank_float_type capital = HANK::ValueNotSet;

		std::vector<hank_float_type> labshareY, labshareN, labfracY, labfracN, labor_occ, wage_occ;

		std::vector<hank_float_type> netwagegrid, yprodgrid;

		int nocc, nprod;

};

class EquilibriumInitial : public Equilibrium {
	public:
		void solve(const Parameters& p, const Model& model);

		void set_from_parameters(const Parameters& p, const Model& model) override;

		void compute_factors(const Model& model) override;
};

class EquilibriumFinal : public Equilibrium {
	public:
		EquilibriumFinal() {}

		EquilibriumFinal(const Equilibrium& other_equm) {
			*this = *(EquilibriumFinal *) &other_equm;
		}

		~EquilibriumFinal() {}

		void solve(const Parameters& p, const Model& model,
			const Equilibrium& initial_equm, const hank_float_type* x);

		void compute_factors(const Model& model) override;
};

class EquilibriumTrans : public Equilibrium {
	public:
		hank_float_type mpshock, pi, pricelev, priceadjust, capadjust, qdot;

		hank_float_type pidot, logydot, elast, firmdiscount, qinvestment, invadjust;

		hank_float_type equity_Adot, equity_Bdot, lIK, bond;

		void set_from_parameters(const Parameters& p, const Model& model) override {
			Equilibrium::set_from_parameters(p, model);
		}

		void compute_factors(const Model& model) override;
};

void solve_trans_equilibrium(std::vector<EquilibriumTrans>& trans_equms,
	const Parameters& p, const Model& model, const EquilibriumInitial& initial_equm,
	const EquilibriumFinal& final_equm, const hank_float_type* deltatransvec);

#endif