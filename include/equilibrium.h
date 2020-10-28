#ifndef _EQUILIBRIUM_H
#define _EQUILIBRIUM_H

#include <hank_config.h>
#include <hank.h>
#include <string>
#include <map>
#include <upwinding.h>

class Parameters;

class Model;

class DistributionStatistics;

class EquilibriumBase : public HankBase
{
	public:
		hank_float_type alpha_Y, alpha_N, price_W, drs_Y, drs_N, riskaver, rho;

		hank_float_type capshareY, capshareN, capfracY, capfracN, capital_Y, capital_N;

		hank_float_type valcapital, labor_Y, labor_N, labor;

		hank_float_type tfp_Y, tfp_N, investment, illprice, illshares, illpricedot;

		hank_float_type netprofit_W, grossprofit_R, netprofit_R, profit;

		hank_float_type rcapital, wage_Y, wage_N, mc_Y, mc_N;

		hank_float_type ra, rb, rnom, pi, dividend_A, dividend_B, equity_A, equity_B;

		hank_float_type Enetwage, taxrev, transfershock, lumptransfer, rborr;

		hank_float_type bond, govbond, labtax, govexp;

		hank_float_type output, varieties, qcapital, priceadjust;

		hank_float_type capital = HANK::ValueNotSet;

		hank_float_type tdelta = 0;

		std::vector<hank_float_type> labshareY, labshareN, labfracY, labfracN, labor_occ, wage_occ;

		std::vector<hank_float_type> netwagegrid;

		hank_float_type depreciation, capadjcost;

		vector3dr V, density;

		Upwinding::Policies policies;

		int nocc, nprod;

		// Variables copied from a Model object
		std::vector<hank_float_type> occYsharegrid, occNsharegrid, occdist, prodgrid, proddist;
};

class Equilibrium : public EquilibriumBase
{
	public:
		virtual std::map<std::string, hank_float_type> variables_map() const override;

		virtual bool is_transition_equilibrium() const {return false;}
};

namespace HANK {
	void print(const Equilibrium& equm);
}

class EquilibriumInitial : public Equilibrium
{
	public:
		void setup(const Parameters& p, const Model& model);

		void solve(const Parameters& p);

		template<typename DistributionStatisticsType, typename HJBType>
		void update_after_solving(const DistributionStatisticsType& stats, const HJBType& hjb);

		std::string title() const override {return "INITIAL EQUILIBRIUM VARIABLES";}
};

class EquilibriumFinal : public Equilibrium
{
	public:
		EquilibriumFinal() {}

		EquilibriumFinal(const EquilibriumBase& other_equm);

		void solve(const Parameters& p, const Equilibrium& initial_equm, const hank_float_type* x);

		std::string title() const override {return "FINAL EQUILIBRIUM VARIABLES";}
};

class EquilibriumTrans : public Equilibrium
{
	public:
		EquilibriumTrans() {}

		EquilibriumTrans(const EquilibriumBase& other_equm);

		hank_float_type mpshock, pricelev, capadjust, qdot;

		hank_float_type pidot, logydot, firmdiscount, qinvestment, invadjust;

		hank_float_type equity_Adot, equity_Bdot, inv_cap_ratio;

		int t = -1;

		std::map<std::string, hank_float_type> variables_map() const override;

		bool is_transition_equilibrium() const override {return true;}

		friend void solve_trans_equilibrium(std::vector<EquilibriumTrans>& trans_equms,
			const Parameters& p, const EquilibriumInitial& initial_equm,
			const EquilibriumFinal& final_equm, const hank_float_type* deltatransvec);

		std::string title() const override {return "TRANSITION EQUILIBRIUM VARIABLES (t = " + std::to_string(t) + ")";}
};

template<typename DistributionStatisticsType, typename HJBType>
void EquilibriumInitial::update_after_solving(const DistributionStatisticsType& stats, const HJBType& hjb)
{
	bond = stats.Eb;
	govbond = equity_B - bond;
	govexp = taxrev + rb * govbond;
	density = stats.density;

	V = hjb.V;
	policies = *hjb.optimal_decisions;
}

#endif