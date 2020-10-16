#ifndef _EQUILIBRIUM_H
#define _EQUILIBRIUM_H

#include <hank_config.h>
#include <hank_types.h>
#include <memory>

class Parameters;

class Model;

template<typename Eltype>
class HankArray {
	public:
		std::unique_ptr<Eltype[]> data = nullptr;

		HankArray() {
			reset(1);
		}

		HankArray(int n) {
			reset(n);
		}

		void reset(int n) {
			data.reset(new Eltype[n]);
			T = n;
		}

		Eltype& operator[](int i) {return data[i];}

		Eltype operator[](int i) const {return data[i];}

		Eltype& get(int i) {return data[i];}

		Eltype get(int i) const {return data[i];}

		int T;
};

class EquilibriumElement {
	public:
		EquilibriumElement() {}

		void create_initial_steady_state(const Parameters& p, const Model& model);

		void set_parameters(const Parameters& p);

		void compute_factors(const Model& model);

		void compute_profits();

		void compute_factor_prices();

		void compute_dividends(const Parameters& p);

		void compute_govt(const Parameters& p, const Model& modelmodel);

		hank_float_type alpha_Y, alpha_N, price_W, drs_Y, drs_N, riskaver, rho;

		hank_float_type capshareY, capshareN, capfracY, capfracN, capital_Y, capital_N;

		hank_float_type labor_Y, labor_N, labor;

		hank_float_type tfp_Y, tfp_N, investment, illprice, illshares, illpricedot;

		hank_float_type netprofit_W, grossprofit_R, netprofit_R, profit;

		hank_float_type rcapital, wage_Y, wage_N, mc_Y, mc_N;

		hank_float_type ra, rb, rnom, pi, dividend_A, dividend_B, equity_A, equity_B;

		hank_float_type Enetwage, taxrev;

		hank_float_type output = 1.0;

		hank_float_type varieties = 1.0;

		hank_float_type qcapital = 1.0;

		hank_float_type capital = HANK::ValueNotSet;

		std::vector<hank_float_type> labshareY, labshareN, labfracY, labfracN, labor_occ, wage_occ;

		std::vector<hank_float_type> netwagegrid, yprodgrid;

		int nocc, nprod;

};

class Equilibrium : public HankArray<EquilibriumElement> {
	public:
		Equilibrium() : HankArray<EquilibriumElement>() {}

		Equilibrium(int n) : HankArray<EquilibriumElement>(n) {}
		// void create_initial_steady_state(const Parameters& p);

		// HankArray<FactorQuantities> compute_factors(const Model& model, const std::vector<hank_float_type>& labor_occ);
};

class TransEquilibriumElement : public EquilibriumElement {
	public:
		hank_float_type mpshock, pi, pricelev, priceadjust;

		hank_float_type pidot, logydot, elast, firmdiscount;
};

class TransEquilibrium : public HankArray<TransEquilibriumElement> {
	public:
		TransEquilibrium() : HankArray<TransEquilibriumElement>() {}

		TransEquilibrium(int n) : HankArray<TransEquilibriumElement>(n) {}
};

#endif