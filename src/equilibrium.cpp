#include <equilibrium.h>
#include <parameters.h>
#include <model.h>
#include <algorithm>
#include <math.h>
#include <distribution_statistics.h>
#include <assert.h>

namespace
{
	double cobb_douglas(hank_float_type capital, hank_float_type labor, hank_float_type alpha);

	void update_equity_variables(const Parameters& p, const std::vector<hank_float_type>& linv, const hank_float_type* deltatransvec,
		const EquilibriumFinal& final_equm, std::vector<EquilibriumTrans>& trans_equms, int T);
}

Equilibrium::Equilibrium() {
	set_pointers();
	HANK::initialize_unset(variable_ptrs);
}

void Equilibrium::set_pointers() {
	variable_ptrs.insert({"output", &output});
	variable_ptrs.insert({"varieties", &varieties});
	variable_ptrs.insert({"qcapital", &qcapital});
	variable_ptrs.insert({"alpha_Y", &alpha_Y});
	variable_ptrs.insert({"alpha_N", &alpha_N});
	variable_ptrs.insert({"price_W", &price_W});
	variable_ptrs.insert({"drs_Y", &drs_Y});
	variable_ptrs.insert({"drs_N", &drs_N});
	variable_ptrs.insert({"riskaver", &riskaver});
	variable_ptrs.insert({"rho", &rho});
	variable_ptrs.insert({"capshareY", &capshareY});
	variable_ptrs.insert({"capshareN", &capshareN});
	variable_ptrs.insert({"capfracY", &capfracY});
	variable_ptrs.insert({"capfracN", &capfracN});
	variable_ptrs.insert({"capital_Y", &capital_Y});
	variable_ptrs.insert({"capital_N", &capital_N});
	variable_ptrs.insert({"valcapital", &valcapital});
	variable_ptrs.insert({"labor_Y", &labor_Y});
	variable_ptrs.insert({"labor_N", &labor_N});
	variable_ptrs.insert({"labor", &labor});
	variable_ptrs.insert({"tfp_Y", &tfp_Y});
	variable_ptrs.insert({"tfp_N", &tfp_N});
	variable_ptrs.insert({"investment", &investment});
	variable_ptrs.insert({"illprice", &illprice});
	variable_ptrs.insert({"illshares", &illshares});
	variable_ptrs.insert({"illpricedot", &illpricedot});
	variable_ptrs.insert({"netprofit_W", &netprofit_W});
	variable_ptrs.insert({"grossprofit_R", &grossprofit_R});
	variable_ptrs.insert({"netprofit_R", &netprofit_R});
	variable_ptrs.insert({"profit", &profit});
	variable_ptrs.insert({"wage_Y", &wage_Y});
	variable_ptrs.insert({"wage_N", &wage_N});
	variable_ptrs.insert({"mc_Y", &mc_Y});
	variable_ptrs.insert({"mc_N", &mc_N});
	variable_ptrs.insert({"ra", &ra});
	variable_ptrs.insert({"rb", &rb});
	variable_ptrs.insert({"rnom", &rnom});
	variable_ptrs.insert({"pi", &pi});
	variable_ptrs.insert({"dividend_A", &dividend_A});
	variable_ptrs.insert({"dividend_B", &dividend_B});
	variable_ptrs.insert({"equity_A", &equity_A});
	variable_ptrs.insert({"equity_B", &equity_B});
	variable_ptrs.insert({"Enetwage", &Enetwage});
	variable_ptrs.insert({"taxrev", &taxrev});
	variable_ptrs.insert({"transfershock", &transfershock});
	variable_ptrs.insert({"lumptransfer", &lumptransfer});
	variable_ptrs.insert({"rborr", &rborr});
	variable_ptrs.insert({"bond", &bond});
	variable_ptrs.insert({"govbond", &govbond});
	variable_ptrs.insert({"labtax", &labtax});
	variable_ptrs.insert({"govexp", &govexp});
	variable_ptrs.insert({"capital", &capital});
}

void Equilibrium::set_from_parameters(const Parameters& p, const Model& model)
{
	alpha_Y = p.alpha_Y;
	alpha_N = p.alpha_N;
	drs_Y = p.drs_Y;
	drs_N = p.drs_N;
	nocc = p.nocc;
	rho = p.rho;
	nprod = model.nprod;
	rborr = p.rborr;
	transfershock = 1.0;
	lumptransfer = p.lumptransfer;
}

void Equilibrium::compute_factors(const Model& model)
{
	// Capital
	capshareY = alpha_Y * price_W * drs_Y;
	capshareN = alpha_N * (1.0 - price_W) * drs_N;
	capfracY = capshareY / (capshareY + capshareN);
	capfracN = capshareN / (capshareY + capshareN);

	// Labor
	labor_Y = 1.0;
	labor_N = 1.0;
	labor = 0;
	labshareY.resize(nocc);
	labshareN.resize(nocc);
	labfracY.resize(nocc);
	labfracN.resize(nocc);
	for (int io=0; io<nocc; ++io) {
		labshareY[io] = (1.0 - alpha_Y) * price_W * drs_Y * model.occYsharegrid[io];
		labshareN[io] = (1.0 - alpha_N) * (1.0 - price_W) * drs_N * model.occNsharegrid[io];
		labfracY[io] = labshareY[io] / (labshareY[io] + labshareN[io]);
		labfracN[io] = labshareN[io] / (labshareY[io] + labshareN[io]);
		labor_Y *= pow(labfracY[io] * labor_occ[io], model.occYsharegrid[io]);
		labor_N *= pow(labfracN[io] * labor_occ[io], model.occNsharegrid[io]);
		labor += labor_occ[io] * model.occdist[io];
	}
}

void Equilibrium::compute_profits()
{
	netprofit_W = price_W * output * (1.0 - drs_Y);
	grossprofit_R = (1.0 - price_W) * output / varieties;
	netprofit_R = varieties * (1.0 - drs_N) * grossprofit_R;
	profit = netprofit_R + netprofit_W;
}

void Equilibrium::compute_factor_prices()
{
	rcapital = (capshareY + capshareN) * output / capital;

	wage_occ.resize(nocc);
	for (int io=0; io<nocc; ++io)
		wage_occ[io] = (labshareN[io] + labshareY[io]) * output / labor_occ[io];

	// Wholesale
	if ( (alpha_Y == 1.0) | (drs_Y == 0.0) )
		wage_Y = 0.0;
	else
		wage_Y = price_W * (1.0 - alpha_Y) * drs_Y * output / labor_Y;

	mc_Y = pow(rcapital / alpha_Y, alpha_Y) * pow(wage_Y / (1.0 - alpha_Y), 1.0 - alpha_Y) / tfp_Y;

	// Expansion
	if ( (alpha_N == 1.0) | (drs_N == 0.0) )
		wage_N = 0.0;
	else
		wage_N = grossprofit_R * (1.0 - alpha_N) * drs_N * varieties / labor_N;

	if ( alpha_N > 0.0 )
		mc_N = pow(rcapital / alpha_N, alpha_N) * pow(wage_N / (1.0 - alpha_N), 1.0 - alpha_N) / tfp_N;
	else
		mc_N = wage_N / tfp_N;
}

void Equilibrium::compute_dividends(const Parameters& p)
{
	ra = rcapital - p.depreciation;
	rb = p.rb;
	dividend_A = p.profdistfracA * profit * (1.0 - p.corptax);
	dividend_B = p.profdistfracB * profit * (1.0 - p.corptax);
	equity_A = dividend_A / ra;
	equity_B = dividend_B / rb;
}

void Equilibrium::compute_netwage(const Parameters& p, const Model& model)
{
	std::vector<hank_float_type> wage_occ_rep;
	Enetwage = 0;

	netwagegrid.resize(nocc * nprod);
	int iy = 0;
	for (int io=0; io<nocc; ++io) {
		for (int ip=0; ip<nprod; ++ip) {
			netwagegrid[iy] = (1.0 - labtax) * model.prodgrid[ip] * wage_occ[io];
			Enetwage += model.occdist[io] * model.proddist[ip] * netwagegrid.back();
			++iy;
		}
	}
}

void EquilibriumInitial::set_from_parameters(const Parameters& p, const Model& model)
{
	Equilibrium::set_from_parameters(p, model);
	price_W = 1.0 - 1.0 / p.elast;
	riskaver = p.riskaver;
	labtax = p.labtax;
	output = 1.0;
	varieties = 1.0;
	qcapital = 1.0;
}

void EquilibriumInitial::solve(const Parameters& p, const Model& model)
{
	set_from_parameters(p, model);

	if ( labor_occ.size() == 0 )
		for (int io=0; io<nocc; ++io)
			labor_occ.push_back(p.hourtarget * p.meanlabeff * model.occdist[io]);

	if ( capital == HANK::ValueNotSet )
		capital = p.target_KY_ratio;

	compute_factors(model);

	tfp_Y = output / pow(cobb_douglas(capital_Y, labor_Y, alpha_Y), drs_Y);
	tfp_N = varieties / pow(cobb_douglas(capital_N, labor_N, alpha_N), drs_N);

	valcapital = capital;
	investment = p.depreciation * capital;

	compute_profits();
	compute_factor_prices();
	compute_dividends(p);
	compute_netwage(p, model);

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

void EquilibriumInitial::compute_factors(const Model& model) {
	Equilibrium::compute_factors(model);
	capital_Y = capfracY * capital;
	capital_N = capfracN * capital;
}

void EquilibriumInitial::update_with_stats(const DistributionStatistics& stats) {
	bond = stats.Eb;
	govbond = equity_B - bond;
	govexp = taxrev + rb * govbond;
}

void EquilibriumInitial::check_results() const {
	HANK::check_if_unset(variable_ptrs);
}

void EquilibriumFinal::set_from_parameters(const Parameters& p, const Model& model)
{
	Equilibrium::set_from_parameters(p, model);
	price_W = 1.0 - 1.0 / p.elast;
	labtax = p.labtax;
	qcapital = 1.0;
}


void EquilibriumFinal::solve(const Parameters& p, const Model& model,
	const Equilibrium& initial_equm, const hank_float_type* x)
{
	set_from_parameters(p, model);

	capital = x[0];

	if ( labor_occ.size() == 0 )
		for (int io=0; io<nocc; ++io)
			labor_occ.push_back(x[io+1]);

	rb = exp(x[nocc+1]);

	compute_factors(model);

	tfp_Y = initial_equm.tfp_Y;
	tfp_N = initial_equm.tfp_N;
	output = tfp_Y * pow(cobb_douglas(capital_Y, labor_Y, alpha_Y), drs_Y);
	varieties = tfp_N * pow(cobb_douglas(capital_N, labor_N, alpha_N), drs_N);

	valcapital = capital;
	investment = p.depreciation * capital;

	compute_profits();
	compute_factor_prices();
	compute_dividends(p);
	compute_netwage(p, model);

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

void EquilibriumFinal::compute_factors(const Model& model) {
	Equilibrium::compute_factors(model);
	capital_Y = capfracY * capital;
	capital_N = capfracN * capital;
}

void EquilibriumFinal::check_results() const {
	HANK::check_if_unset(variable_ptrs);
}

EquilibriumTrans::EquilibriumTrans() : Equilibrium() {
	set_pointers();
	HANK::initialize_unset(variable_ptrs);
}

EquilibriumTrans::EquilibriumTrans(const Equilibrium& other_equm) : Equilibrium() {
	set_pointers();
	auto variable_ptrs_copy = variable_ptrs;

	HANK::initialize_unset(variable_ptrs);
	*this = *(EquilibriumTrans *) &other_equm;
	variable_ptrs = variable_ptrs_copy;

	labshareY.clear();
	labshareN.clear();
	labfracY.clear();
	labfracN.clear();
	labor_occ.clear();
	wage_occ.clear();
	netwagegrid.clear();
	yprodgrid.clear();
}

void EquilibriumTrans::set_pointers() {
	variable_ptrs.insert({"mpshock", &mpshock});
	variable_ptrs.insert({"pricelev", &pricelev});
	variable_ptrs.insert({"priceadjust", &priceadjust});
	variable_ptrs.insert({"capadjust", &capadjust});
	variable_ptrs.insert({"qdot", &qdot});
	variable_ptrs.insert({"pidot", &pidot});
	variable_ptrs.insert({"logydot", &logydot});
	variable_ptrs.insert({"firmdiscount", &firmdiscount});
	variable_ptrs.insert({"qinvestment", &qinvestment});
	variable_ptrs.insert({"invadjust", &invadjust});
	variable_ptrs.insert({"equity_Adot", &equity_Adot});
	variable_ptrs.insert({"equity_Bdot", &equity_Bdot});
	variable_ptrs.insert({"inv_cap_ratio", &inv_cap_ratio});
	variable_ptrs.insert({"output", &output});
}

void EquilibriumTrans::compute_factors(const Model& model) {
	Equilibrium::compute_factors(model);
	capital_Y = pow(output / tfp_Y, 1.0 / drs_Y) / pow(labor_Y, 1.0 - alpha_Y);
	capital_Y = pow(capital_Y, 1.0 / alpha_Y);
	capital = capital_Y / capfracY;
	capital_N = capfracN * capital;
}

void EquilibriumTrans::check_results() const {
	HANK::check_if_unset(variable_ptrs);
}

void solve_trans_equilibrium(std::vector<EquilibriumTrans>& trans_equms,
	const Parameters& p, const Model& model, const EquilibriumInitial& initial_equm,
	const EquilibriumFinal& final_equm, const hank_float_type* deltatransvec)
{
	int T = trans_equms.size();
	std::vector<hank_float_type> linv(T);

	for (int it=0; it<T; ++it) {
		double deltatrans = deltatransvec[it];
		trans_equms[it].compute_factors(model);

		trans_equms[it].varieties = trans_equms[it].tfp_N * pow(
			cobb_douglas(trans_equms[it].capital_N, trans_equms[it].labor_N, trans_equms[it].alpha_N),
			trans_equms[it].drs_N);

		trans_equms[it].compute_profits();
		trans_equms[it].compute_factor_prices();
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

			trans_equms[it].inv_cap_ratio = model.capadjcost1inv(trans_equms[it].qcapital - 1.0);

			if (it < T - 1)
				trans_equms[it].qdot = (trans_equms[it+1].qcapital - trans_equms[it].qcapital) / deltatrans;
			else
				trans_equms[it].qdot = 0;

			trans_equms[it].capadjust = model.capadjcost(trans_equms[it].inv_cap_ratio);
			trans_equms[it].ra = (trans_equms[it].rcapital + trans_equms[it].inv_cap_ratio * model.capadjcost1(trans_equms[it].inv_cap_ratio) - trans_equms[it].capadjust + trans_equms[it].qdot)
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

		trans_equms[T-1].inv_cap_ratio = model.capadjcost1inv(trans_equms[T-1].qcapital - 1.0);
		for (int it=0; it<T; ++it)
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
		trans_equms[it].compute_netwage(p, model);
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