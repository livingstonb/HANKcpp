#ifndef _EQUILIBRIUM_H
#define _EQUILIBRIUM_H

#include <hank_config.h>
#include <hank.h>
#include <string>
#include <map>

class Parameters;

class Model;

class DistributionStatistics;

class EquilibriumBase
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

		hank_float_type output, varieties, qcapital;

		hank_float_type capital = HANK::ValueNotSet;

		std::vector<hank_float_type> labshareY, labshareN, labfracY, labfracN, labor_occ, wage_occ;

		std::vector<hank_float_type> netwagegrid;

		hank_float_type depreciation, capadjcost;

		int nocc, nprod;

		// Variables copied from a Model object
		std::vector<hank_float_type> occYsharegrid, occNsharegrid, occdist, prodgrid, proddist;
};

class Equilibrium : public EquilibriumBase
{
	public:
		virtual std::map<std::string, hank_float_type> variables_map() const;
};

namespace HANK {
	void print(const Equilibrium& equm);
}

class EquilibriumInitial : public Equilibrium
{
	public:
		void setup(const Parameters& p, const Model& model);

		void solve(const Parameters& p);

		template<typename DistributionStatisticsType>
		void update_with_stats(const DistributionStatisticsType& stats);
};

class EquilibriumFinal : public Equilibrium
{
	public:
		EquilibriumFinal() {}

		EquilibriumFinal(const EquilibriumBase& other_equm);

		void solve(const Parameters& p, const Equilibrium& initial_equm, const hank_float_type* x);
};

class EquilibriumTrans : public Equilibrium
{
	public:
		EquilibriumTrans() {}

		EquilibriumTrans(const EquilibriumBase& other_equm);

		hank_float_type mpshock, pricelev, priceadjust, capadjust, qdot;

		hank_float_type pidot, logydot, firmdiscount, qinvestment, invadjust;

		hank_float_type equity_Adot, equity_Bdot, inv_cap_ratio;

		std::map<std::string, hank_float_type> variables_map() const override;

		friend void solve_trans_equilibrium(std::vector<EquilibriumTrans>& trans_equms,
			const Parameters& p, const EquilibriumInitial& initial_equm,
			const EquilibriumFinal& final_equm, const hank_float_type* deltatransvec);
};

template<typename DistributionStatisticsType>
void EquilibriumInitial::update_with_stats(const DistributionStatisticsType& stats)
{
	bond = stats.Eb;
	govbond = equity_B - bond;
	govexp = taxrev + rb * govbond;
}

namespace HANK {
	void print(const Equilibrium& equm);
}

#endif