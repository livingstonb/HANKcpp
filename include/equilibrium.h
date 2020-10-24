#ifndef _EQUILIBRIUM_H
#define _EQUILIBRIUM_H

#include <hank_config.h>
#include <hank_types.h>
#include <string>
#include <map>

class Parameters;

class DistributionStatistics;

class Equilibrium {
	private:
		void set_from_parameters(const Parameters& p);

		template<typename ModelType>
		void set_from_model(const ModelType& model);

	public:
		Equilibrium() {};

		void create_transition();

		template<typename ModelType>
		void setup(const Parameters& p, const ModelType& model);

		virtual void compute_factors();

		virtual std::map<std::string, hank_float_type> get_variables_map() const;

		virtual void print() const;

		void compute_profits();

		void compute_factor_prices();

		void compute_dividends(const Parameters& p);

		void compute_netwage(const Parameters& p);

		hank_float_type alpha_Y, alpha_N, price_W, drs_Y, drs_N, riskaver, rho;

		hank_float_type capshareY, capshareN, capfracY, capfracN, capital_Y, capital_N;

		hank_float_type valcapital, labor_Y, labor_N, labor;

		hank_float_type tfp_Y, tfp_N, investment, illprice, illshares, illpricedot;

		hank_float_type netprofit_W, grossprofit_R, netprofit_R, profit;

		hank_float_type rcapital, wage_Y, wage_N, mc_Y, mc_N;

		hank_float_type ra, rb, rnom, pi, dividend_A, dividend_B, equity_A, equity_B;

		hank_float_type Enetwage, taxrev, transfershock, lumptransfer, rborr;

		hank_float_type bond, govbond, labtax, govexp;

		hank_float_type output, varieties, qcapital;

		hank_float_type capital = HANK::ValueNotSet;

		std::vector<hank_float_type> labshareY, labshareN, labfracY, labfracN, labor_occ, wage_occ;

		std::vector<hank_float_type> netwagegrid;

		hank_float_type depreciation, capadjcost;

		int nocc, nprod;

		// Variables copied from a Model object
		std::vector<hank_float_type> occYsharegrid, occNsharegrid, occdist, prodgrid, proddist;
};

class EquilibriumInitial : public Equilibrium {
	public:
		void solve(const Parameters& p);

		void compute_factors() override;

		template<typename T>
		void update_with_stats(const T& stats);

		void print() const override;
};

class EquilibriumFinal : public Equilibrium {
	public:
		EquilibriumFinal() {}

		EquilibriumFinal(const Equilibrium& other_equm);

		void solve(const Parameters& p, const Equilibrium& initial_equm, const hank_float_type* x);

		void compute_factors() override;

		void print() const override;
};

class EquilibriumTrans : public Equilibrium {
	public:
		EquilibriumTrans() {};

		EquilibriumTrans(const Equilibrium& other_equm);

		hank_float_type mpshock, pricelev, priceadjust, capadjust, qdot;

		hank_float_type pidot, logydot, firmdiscount, qinvestment, invadjust;

		hank_float_type equity_Adot, equity_Bdot, inv_cap_ratio;

		void compute_factors() override;

		void print() const override;

		virtual std::map<std::string, hank_float_type> get_variables_map() const override;
};

template<typename ModelType>
void Equilibrium::setup(const Parameters& p, const ModelType& model)
{
	set_from_parameters(p);
	set_from_model(model);
}

template<typename ModelType>
void Equilibrium::set_from_model(const ModelType& model)
{
	nprod = model.nprod;
	occYsharegrid = model.occYsharegrid;
	occNsharegrid = model.occNsharegrid;
	occdist = model.occdist;
	prodgrid = model.prodgrid;
	proddist = model.proddist;
}

template<typename DistributionStatisticsType>
void EquilibriumInitial::update_with_stats(const DistributionStatisticsType& stats)
{
	bond = stats.Eb;
	govbond = equity_B - bond;
	govexp = taxrev + rb * govbond;
}

void solve_trans_equilibrium(std::vector<EquilibriumTrans>& trans_equms,
	const Parameters& p, const EquilibriumInitial& initial_equm,
	const EquilibriumFinal& final_equm, const hank_float_type* deltatransvec);

#endif