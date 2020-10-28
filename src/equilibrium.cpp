#include <equilibrium.h>
#include <parameters.h>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <model_functions.h>
#include <model.h>

namespace
{
	double cobb_douglas(hank_float_type capital, hank_float_type labor, hank_float_type alpha);

	void update_equity_variables(const Parameters& p, const std::vector<hank_float_type>& linv, const hank_float_type* deltatransvec,
		const EquilibriumFinal& final_equm, std::vector<EquilibriumTrans>& trans_equms, int T);

	void set_from_parameters(EquilibriumInitial* equm, const Parameters& p);

	void set_from_model(EquilibriumInitial* equm, const Model& model);

	void compute_profits(Equilibrium* equm);

	void compute_factor_prices(Equilibrium* equm);

	void compute_dividends(Equilibrium* equm, const Parameters& p);

	void compute_netwage(Equilibrium* equm);

	void compute_factors(Equilibrium* equm);

	void compute_factors(EquilibriumInitial* equm);

	void compute_factors(EquilibriumFinal* equm);

	void compute_factors(EquilibriumTrans* equm);
}

namespace HANK {
	void print(const Equilibrium& equm)
	{
		std::map<std::string, hank_float_type> variables = equm.variables_map();
		print(variables, "EQUILIBRIUM VARIABLES");
	}
}

std::map<std::string, hank_float_type> Equilibrium::variables_map() const
{
	std::map<std::string, hank_float_type> variables;

	variables.insert({"alpha_Y", alpha_Y});
	variables.insert({"alpha_N", alpha_N});
	variables.insert({"drs_Y", drs_Y});
	variables.insert({"drs_N", drs_N});
	variables.insert({"rb", rb});
	variables.insert({"rborr", rborr});
	variables.insert({"transfershock", transfershock});
	variables.insert({"lumptransfer", lumptransfer});
	variables.insert({"price_W", price_W});
	variables.insert({"riskaver", riskaver});
	variables.insert({"rho", rho});
	variables.insert({"labtax", labtax});
	variables.insert({"output", output});
	variables.insert({"varieties", varieties});
	variables.insert({"qcapital", qcapital});
	variables.insert({"capshareY", capshareY});
	variables.insert({"capshareN", capshareN});
	variables.insert({"capfracY", capfracY});
	variables.insert({"capfracN", capfracN});
	variables.insert({"capital_Y", capital_Y});
	variables.insert({"capital_N", capital_N});
	variables.insert({"valcapital", valcapital});
	variables.insert({"labor_Y", labor_Y});
	variables.insert({"labor_N", labor_N});
	variables.insert({"labor", labor});
	variables.insert({"tfp_Y", tfp_Y});
	variables.insert({"tfp_N", tfp_N});
	variables.insert({"investment", investment});
	variables.insert({"illprice", illprice});
	variables.insert({"illshares", illshares});
	variables.insert({"illpricedot", illpricedot});
	variables.insert({"netprofit_W", netprofit_W});
	variables.insert({"grossprofit_R", grossprofit_R});
	variables.insert({"netprofit_R", netprofit_R});
	variables.insert({"profit", profit});
	variables.insert({"rcapital", rcapital});
	variables.insert({"wage_Y", wage_Y});
	variables.insert({"wage_N", wage_N});
	variables.insert({"mc_Y", mc_Y});
	variables.insert({"mc_N", mc_N});
	variables.insert({"ra", ra});
	variables.insert({"rnom", rnom});
	variables.insert({"pi", pi});
	variables.insert({"dividend_A", dividend_A});
	variables.insert({"dividend_B", dividend_B});
	variables.insert({"equity_A", equity_A});
	variables.insert({"equity_B", equity_B});
	variables.insert({"Enetwage", Enetwage});
	variables.insert({"taxrev", taxrev});
	variables.insert({"bond", bond});
	variables.insert({"govbond", govbond});
	variables.insert({"govexp", govexp});
	variables.insert({"capital", capital});
	variables.insert({"priceadjust", priceadjust});

	return variables;
}

void EquilibriumInitial::setup(const Parameters& p, const Model& model)
{
	priceadjust = 0;
	set_from_parameters(this, p);
	set_from_model(this, model);
}

void EquilibriumInitial::solve(const Parameters& p)
{
	if ( labor_occ.size() == 0 )
		for (int io=0; io<nocc; ++io)
			labor_occ.push_back(p.hourtarget * p.meanlabeff * occdist[io]);

	if ( capital == HANK::ValueNotSet )
		capital = p.target_KY_ratio;

	compute_factors(this);

	tfp_Y = output / pow(cobb_douglas(capital_Y, labor_Y, alpha_Y), drs_Y);
	tfp_N = varieties / pow(cobb_douglas(capital_N, labor_N, alpha_N), drs_N);

	valcapital = capital;
	investment = depreciation * capital;

	compute_profits(this);
	compute_factor_prices(this);
	compute_dividends(this, p);
	compute_netwage(this);

	taxrev = p.corptax * profit - p.lumptransfer;
	for (int io=0; io<nocc; ++io)
		taxrev += labtax * wage_occ[io] * labor_occ[io];

	if ( p.taxHHProfitIncome )
		taxrev += labtax * p.profdistfracW * profit * (1.0 - p.corptax);

	pi = p.pi;
	rnom = rb - p.pi;
	illprice = 1;
	illpricedot = 0;
	illshares = capital + equity_A;
}

EquilibriumFinal::EquilibriumFinal(const EquilibriumBase& other_equm)
{
	EquilibriumBase* eqbase = (EquilibriumBase *) this;
	*eqbase = other_equm;

	// labshareY.clear();
	// labshareN.clear();
	// labfracY.clear();
	// labfracN.clear();
	// labor_occ.clear();
	// wage_occ.clear();
	// netwagegrid.clear();
}

void EquilibriumFinal::solve(const Parameters& p, const Equilibrium& initial_equm, const hank_float_type* x)
{
	capital = x[0];

	if ( labor_occ.size() == 0 )
		for (int io=0; io<nocc; ++io)
			labor_occ.push_back(x[io+1]);

	rb = exp(x[nocc+1]);

	compute_factors(this);

	tfp_Y = initial_equm.tfp_Y;
	tfp_N = initial_equm.tfp_N;
	output = tfp_Y * pow(cobb_douglas(capital_Y, labor_Y, alpha_Y), drs_Y);
	varieties = tfp_N * pow(cobb_douglas(capital_N, labor_N, alpha_N), drs_N);

	valcapital = capital;
	investment = depreciation * capital;

	compute_profits(this);
	compute_factor_prices(this);
	compute_dividends(this, p);
	compute_netwage(this);

	taxrev = p.corptax * profit - p.lumptransfer;
	for (int io=0; io<nocc; ++io)
		taxrev += labtax * wage_occ[io] * labor_occ[io];

	if ( p.taxHHProfitIncome )
		taxrev += labtax * p.profdistfracW * profit * (1.0 - p.corptax);

	pi = p.pi;
	rnom = rb - p.pi;
	illprice = capital + equity_A;
	illpricedot = 0;
	illshares = initial_equm.illshares;

	govbond = (govexp - taxrev) / rb;
	bond = equity_B - govbond;
}

EquilibriumTrans::EquilibriumTrans(const EquilibriumBase& other_equm)
{
	EquilibriumBase* eqbase = (EquilibriumBase *) this;
	*eqbase = other_equm;

	// labshareY.clear();
	// labshareN.clear();
	// labfracY.clear();
	// labfracN.clear();
	// wage_occ.clear();
	// netwagegrid.clear();
}

std::map<std::string, hank_float_type> EquilibriumTrans::variables_map() const
{
	std::map<std::string, hank_float_type> variables = Equilibrium::variables_map();

	variables.insert({"mpshock", mpshock});
	variables.insert({"pricelev", pricelev});
	variables.insert({"capadjust", capadjust});
	variables.insert({"qdot", qdot});
	variables.insert({"pidot", pidot});
	variables.insert({"logydot", logydot});
	variables.insert({"firmdiscount", firmdiscount});
	variables.insert({"qinvestment", qinvestment});
	variables.insert({"invadjust", invadjust});
	variables.insert({"equity_Adot", equity_Adot});
	variables.insert({"equity_Bdot", equity_Bdot});

	return variables;
}

void solve_trans_equilibrium(std::vector<EquilibriumTrans>& trans_equms,
	const Parameters& p, const EquilibriumInitial& initial_equm,
	const EquilibriumFinal& final_equm, const hank_float_type* deltatransvec)
{
	int T = trans_equms.size();
	std::vector<hank_float_type> linv(T);

	for (int it=0; it<T; ++it) {
		double deltatrans = deltatransvec[it];
		compute_factors(&trans_equms[it]);

		trans_equms[it].varieties = trans_equms[it].tfp_N * pow(
			cobb_douglas(trans_equms[it].capital_N, trans_equms[it].labor_N, trans_equms[it].alpha_N),
			trans_equms[it].drs_N);

		compute_profits(&trans_equms[it]);
		compute_factor_prices(&trans_equms[it]);
		trans_equms[it].investment = p.depreciation * trans_equms[it].capital;
		if (it < T - 1)
			trans_equms[it].investment += (trans_equms[it+1].capital - trans_equms[it].capital) / deltatrans;


		// Value of capital and ra
		if ( (p.capadjcost == 0) & (p.invadjcost == 0) ) {
			trans_equms[it].capadjust = 0;
			trans_equms[it].qdot = 0;
			trans_equms[it].ra = trans_equms[it].rcapital - p.depreciation;
			trans_equms[it].qinvestment = 0;
			linv[it] = 0;
			trans_equms[it].inv_cap_ratio = 0;
		}
		else if ( p.capadjcost > 0 ) {
			trans_equms[it].qinvestment = 0;
			trans_equms[it].invadjust = 0;

			trans_equms[it].inv_cap_ratio = ModelFunctions::capadjcost1inv(
				trans_equms[it].qcapital - 1.0, p.capadjcost, p.depreciation);

			if (it < T - 1)
				trans_equms[it].qdot = (trans_equms[it+1].qcapital - trans_equms[it].qcapital) / deltatrans;
			else
				trans_equms[it].qdot = 0;

			trans_equms[it].capadjust = ModelFunctions::capadjcost1(trans_equms[it].inv_cap_ratio, p.capadjcost, p.depreciation);
			trans_equms[it].ra =
				(trans_equms[it].rcapital + trans_equms[it].inv_cap_ratio * trans_equms[it].capadjust - trans_equms[it].capadjust + trans_equms[it].qdot)
				/ trans_equms[it].qcapital - p.depreciation;

			linv[it] = 0;
		}
	}

	
	if ( p.invadjcost > 0 ) {
		trans_equms[T-1].qdot = 0;
		trans_equms[T-1].qinvestment = 0;
		trans_equms[T-1].ra = (trans_equms[T-1].rcapital + trans_equms[T-1].qdot) / trans_equms[T-1].qcapital - p.depreciation;

		for (int it=T-2; it>=0; --it) {
			trans_equms[it].qdot = (trans_equms[it+1].qcapital - trans_equms[it].qcapital) / deltatransvec[it];
			trans_equms[it].ra = (trans_equms[it].rcapital + trans_equms[it].qdot) / trans_equms[it].qcapital - p.depreciation;

			double qinv1 = trans_equms[it].qcapital - 1.0 - pow(trans_equms[it+1].qinvestment, 2) / (2.0 * p.invadjcost);
			double qinv2 = 1.0 + (trans_equms[it].ra - trans_equms[it+1].qinvestment / p.invadjcost) * deltatransvec[it];
			trans_equms[it].qinvestment = (qinv1 * deltatransvec[it] + trans_equms[it+1].qinvestment) / qinv2;
		}

		linv[0] = initial_equm.investment;
		for (int it=1; it<T; ++it)
			linv[it] = linv[it-1] / (1.0 - deltatransvec[it] * trans_equms[it-1].qinvestment / p.invadjcost);

		trans_equms[T-1].inv_cap_ratio = ModelFunctions::capadjcost1inv(trans_equms[T-1].qcapital - 1.0, p.capadjcost, p.depreciation);
		for (int it=0; it<T-1; ++it)
			trans_equms[it].inv_cap_ratio = linv[it] / trans_equms[it].capital;
	}

	update_equity_variables(p, linv, deltatransvec, final_equm, trans_equms, T);

	// Government
	std::vector<hank_float_type> grosslabtaxinc(T);
	for (int it=0; it<T; ++it) {
		grosslabtaxinc[it] = 0;
		for (int io=0; io<p.nocc; ++io)
			grosslabtaxinc[it] += trans_equms[it].wage_occ[io] * trans_equms[it].labor_occ[io];
	}

	if (p.adjGovBudgetConstraint == GovBCAdjType::fiscal) {
		trans_equms[0].govbond = initial_equm.govbond;

		for (int it=0; it<T; ++it) {
			trans_equms[it].labtax = initial_equm.labtax;
			trans_equms[it].govexp = initial_equm.govexp;
			double transfershock = initial_equm.lumptransfer * trans_equms[it].transfershock - initial_equm.lumptransfer;
			double grosstaxrev = trans_equms[it].labtax * grosslabtaxinc[it] + p.corptax * trans_equms[it].profit - transfershock;

			if ( p.taxHHProfitIncome )
				grosstaxrev += trans_equms[it].labtax * p.profdistfracW * trans_equms[it].profit * (1.0 - p.corptax);
		

			// Solve gov budget constraint forwards
			if ( it < T - 1 ) {
				trans_equms[it+1].govbond = (trans_equms[it].govbond + deltatransvec[it]
					* ((trans_equms[it].rb + p.govdebtadjwedge) * initial_equm.govbond + grosstaxrev - trans_equms[it].govexp - initial_equm.lumptransfer))
					/ (1.0 + deltatransvec[it] * p.govdebtadjwedge);
			}

			trans_equms[it].lumptransfer = initial_equm.lumptransfer
				+ (trans_equms[it].rb + p.govdebtadjwedge)  * (trans_equms[it].govbond - initial_equm.govbond) + transfershock;
			trans_equms[it].taxrev = grosstaxrev + transfershock - trans_equms[it].lumptransfer;

			if ( trans_equms[it].taxrev < 0 )
				std::cout << "Warning: negative transfers\n";
		}
	}

	for (int it=0; it<T; ++it) {
		compute_netwage(&trans_equms[it]);
		trans_equms[it].bond = trans_equms[it].equity_B - trans_equms[it].govbond;
		trans_equms[it].rborr = trans_equms[it].rb + p.borrwedge;
	}
}

namespace
{
	double cobb_douglas(hank_float_type capital, hank_float_type labor, hank_float_type alpha)
	{
		return pow(capital, alpha) * pow(labor, 1.0 - alpha);
	}

	void set_from_parameters(EquilibriumInitial* equm, const Parameters& p)
	{
		equm->alpha_Y = p.alpha_Y;
		equm->alpha_N = p.alpha_N;
		equm->drs_Y = p.drs_Y;
		equm->drs_N = p.drs_N;
		equm->nocc = p.nocc;
		equm->rho = p.rho;
		equm->rb = p.rb;
		equm->rborr = p.rborr;
		equm->transfershock = 1.0;
		equm->lumptransfer = p.lumptransfer;
		equm->price_W = 1.0 - 1.0 / p.elast;
		equm->riskaver = p.riskaver;
		equm->labtax = p.labtax;
		equm->output = 1.0;
		equm->varieties = 1.0;
		equm->qcapital = 1.0;

		equm->capadjcost = p.capadjcost;
		equm->depreciation = p.depreciation;
	}

	void set_from_model(EquilibriumInitial* equm, const Model& model)
	{
		equm->nprod = model.nprod;
		equm->occYsharegrid = model.occYsharegrid;
		equm->occNsharegrid = model.occNsharegrid;
		equm->occdist = model.occdist;
		equm->prodgrid = model.prodgrid;
		equm->proddist = model.proddist;
	}

	void compute_profits(Equilibrium* equm)
	{
		equm->netprofit_W = equm->price_W * equm->output * (1.0 - equm->drs_Y);
		equm->grossprofit_R = (1.0 - equm->price_W) * equm->output / equm->varieties;
		equm->netprofit_R = equm->varieties * (1.0 - equm->drs_N) * equm->grossprofit_R;
		equm->profit = equm->netprofit_R + equm->netprofit_W - equm->priceadjust;
	}

	void compute_factor_prices(Equilibrium* equm)
	{
		equm->rcapital = (equm->capshareY + equm->capshareN) * equm->output / equm->capital;

		equm->wage_occ.resize(equm->nocc);
		for (int io=0; io<equm->nocc; ++io)
			equm->wage_occ[io] = (equm->labshareN[io] + equm->labshareY[io]) * equm->output / equm->labor_occ[io];

		// Wholesale
		if ( (equm->alpha_Y == 1.0) | (equm->drs_Y == 0.0) )
			equm->wage_Y = 0.0;
		else
			equm->wage_Y = equm->price_W * (1.0 - equm->alpha_Y) * equm->drs_Y * equm->output / equm->labor_Y;

		equm->mc_Y = pow(equm->rcapital / equm->alpha_Y, equm->alpha_Y)
			* pow(equm->wage_Y / (1.0 - equm->alpha_Y), 1.0 - equm->alpha_Y) / equm->tfp_Y;

		// Expansion
		if ( (equm->alpha_N == 1.0) | (equm->drs_N == 0.0) )
			equm->wage_N = 0.0;
		else
			equm->wage_N = equm->grossprofit_R * (1.0 - equm->alpha_N) * equm->drs_N * equm->varieties / equm->labor_N;

		if ( equm->alpha_N > 0.0 )
			equm->mc_N = pow(equm->rcapital / equm->alpha_N, equm->alpha_N)
			* pow(equm->wage_N / (1.0 - equm->alpha_N), 1.0 - equm->alpha_N) / equm->tfp_N;
		else
			equm->mc_N = equm->wage_N / equm->tfp_N;
	}

	void compute_dividends(Equilibrium* equm, const Parameters& p)
	{
		equm->ra = equm->rcapital - equm->depreciation;
		equm->dividend_A = p.profdistfracA * equm->profit * (1.0 - p.corptax);
		equm->dividend_B = p.profdistfracB * equm->profit * (1.0 - p.corptax);
		equm->equity_A = equm->dividend_A / equm->ra;
		equm->equity_B = equm->dividend_B / equm->rb;
	}

	void compute_netwage(Equilibrium* equm)
	{
		equm->Enetwage = 0;
		equm->netwagegrid.resize(equm->nocc * equm->nprod);
		int iy = 0;
		for (int io=0; io<equm->nocc; ++io) {
			for (int ip=0; ip<equm->nprod; ++ip) {
				equm->netwagegrid[iy] = (1.0 - equm->labtax) * equm->prodgrid[ip] * equm->wage_occ[io];
				equm->Enetwage += equm->occdist[io] * equm->proddist[ip] * equm->netwagegrid[iy];
				++iy;
			}
		}
	}

	void compute_factors(Equilibrium* equm)
	{
		// Capital
		equm->capshareY = equm->alpha_Y * equm->price_W * equm->drs_Y;
		equm->capshareN = equm->alpha_N * (1.0 - equm->price_W) * equm->drs_N;
		equm->capfracY = equm->capshareY / (equm->capshareY + equm->capshareN);
		equm->capfracN = equm->capshareN / (equm->capshareY + equm->capshareN);

		// Labor
		equm->labor_Y = 1.0;
		equm->labor_N = 1.0;
		equm->labor = 0;
		equm->labshareY.resize(equm->nocc);
		equm->labshareN.resize(equm->nocc);
		equm->labfracY.resize(equm->nocc);
		equm->labfracN.resize(equm->nocc);
		for (int io=0; io<equm->nocc; ++io) {
			equm->labshareY[io] = (1.0 - equm->alpha_Y) * equm->price_W * equm->drs_Y * equm->occYsharegrid[io];
			equm->labshareN[io] = (1.0 - equm->alpha_N) * (1.0 - equm->price_W) * equm->drs_N * equm->occNsharegrid[io];
			equm->labfracY[io] = equm->labshareY[io] / (equm->labshareY[io] + equm->labshareN[io]);
			equm->labfracN[io] = equm->labshareN[io] / (equm->labshareY[io] + equm->labshareN[io]);
			equm->labor_Y *= pow(equm->labfracY[io] * equm->labor_occ[io], equm->occYsharegrid[io]);
			equm->labor_N *= pow(equm->labfracN[io] * equm->labor_occ[io], equm->occNsharegrid[io]);
			equm->labor += equm->labor_occ[io] * equm->occdist[io];
		}
	}

	void compute_factors(EquilibriumInitial* equm)
	{
		compute_factors((Equilibrium *) equm);
		equm->capital_Y = equm->capfracY * equm->capital;
		equm->capital_N = equm->capfracN * equm->capital;
	}

	void compute_factors(EquilibriumFinal* equm)
	{
		compute_factors((Equilibrium *) equm);
		equm->capital_Y = equm->capfracY * equm->capital;
		equm->capital_N = equm->capfracN * equm->capital;
	}

	void compute_factors(EquilibriumTrans* equm)
	{
		compute_factors((Equilibrium *) equm);
		equm->capital_Y = pow(equm->output / equm->tfp_Y, 1.0 / equm->drs_Y) / pow(equm->labor_Y, 1.0 - equm->alpha_Y);
		equm->capital_Y = pow(equm->capital_Y, 1.0 / equm->alpha_Y);
		equm->capital = equm->capital_Y / equm->capfracY;
		equm->capital_N = equm->capfracN * equm->capital;
	}

	void update_equity_variables(const Parameters& p, const std::vector<hank_float_type>& linv, const hank_float_type* deltatransvec,
		const EquilibriumFinal& final_equm, std::vector<EquilibriumTrans>& trans_equms, int T)
	{
		for (int it=0; it<T; ++it) {
			trans_equms[it].valcapital = trans_equms[it].qcapital * trans_equms[it].capital + trans_equms[it].qinvestment * linv[it];
			trans_equms[it].dividend_A = p.profdistfracA * trans_equms[it].profit * (1.0 - p.corptax);
			trans_equms[it].dividend_B = p.profdistfracB * trans_equms[it].profit * (1.0 - p.corptax);
		}

		trans_equms[T-1].equity_A = final_equm.equity_A;
		trans_equms[T-1].equity_B = final_equm.equity_B;
		trans_equms[T-1].equity_Adot = 0;
		trans_equms[T-1].equity_Bdot = 0;

		for (int it=T-2; it>=0; --it) {
			trans_equms[it].equity_A = (trans_equms[it+1].equity_A + trans_equms[it].dividend_A * deltatransvec[it])
				/ (1.0 + deltatransvec[it] * trans_equms[it].ra);
			trans_equms[it].equity_B = (trans_equms[it+1].equity_B + trans_equms[it].dividend_B * deltatransvec[it])
				/ (1.0 + deltatransvec[it] * trans_equms[it].rb);
			trans_equms[it].equity_Adot = (trans_equms[it+1].equity_A - trans_equms[it].equity_A) / deltatransvec[it];
			trans_equms[it].equity_Bdot = (trans_equms[it+1].equity_B - trans_equms[it].equity_B) / deltatransvec[it];
		}

		for (int it=0; it<T; ++it) {
			trans_equms[it].illprice = (trans_equms[it].valcapital + trans_equms[it].equity_A) / trans_equms[it].illshares;
		}

		trans_equms[T-1].illpricedot = 0;
		for (int it=T-2; it>=0; --it)
			trans_equms[it].illpricedot = (trans_equms[it+1].illprice - trans_equms[it].illprice) / deltatransvec[it];
	}
}