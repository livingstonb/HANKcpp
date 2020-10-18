#include <equilibrium.h>
#include <parameters.h>
#include <model.h>
#include <algorithm>

namespace
{
	double cobb_douglas(hank_float_type capital, hank_float_type labor, hank_float_type alpha);

	void update_equity_variables(const Parameters& p, const std::vector<hank_float_type>& linv, const hank_float_type* deltatransvec,
		const EquilibriumFinal& final_equm, std::vector<EquilibriumTrans>& trans_equms, int T);
}

void Equilibrium::set_from_parameters(const Parameters& p, const Model& model)
{
	alpha_Y = p.alpha_Y;
	alpha_N = p.alpha_N;
	drs_Y = p.drs_Y;
	drs_N = p.drs_N;
	nocc = p.nocc;
	rho = p.rho;

	if ( is_initial_steady_state() ) {
		price_W = 1.0 - 1.0 / p.elast;
		riskaver = p.riskaver;
	}

	nprod = model.nprod;
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
	for (int io=0; io<nocc; ++io) {
		labshareY.push_back((1.0 - alpha_Y) * price_W * drs_Y * model.occYsharegrid[io]);
		labshareN.push_back((1.0 - alpha_N) * (1.0 - price_W) * drs_N * model.occNsharegrid[io]);
		labfracY.push_back(labshareY[io] / (labshareY[io] + labshareN[io]));
		labfracN.push_back(labshareN[io] / (labshareY[io] + labshareN[io]));
		labor_Y *= pow(labfracY[io] * labor_occ[io], model.occYsharegrid[io]);
		labor_N *= pow(labfracN[io] * labor_occ[io], model.occNsharegrid[io]);
		labor += labor_occ[io] * model.occdist[io];
	}

	if ( is_initial_steady_state() | is_final_steady_state() ) {
		capital_Y = capfracY * capital;
		capital_N = capfracN * capital;
	}
	else if ( is_trans_equilibrium() ) {
		capital_Y = pow(output / tfp_Y, 1.0 / drs_Y) / pow(labor_Y, 1.0 - alpha_Y);
		capital_Y = pow(capital_Y, 1.0 / alpha_Y);
		capital = capital_Y / capfracY;
		capital_N = capfracN * capital;
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

	for (int io=0; io<nocc; ++io)
		wage_occ.push_back((labshareN[io] + labshareY[io]) * output / labor_occ[io]);

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

void Equilibrium::compute_govt(const Parameters& p, const Model& model)
{
	std::vector<hank_float_type> wage_occ_rep;
	Enetwage = 0;

	for (int io=0; io<nocc; ++io) {
		for (int ip=0; ip<model.nprod; ++ip) {
			netwagegrid.push_back((1.0 - p.labtax) * model.prodgrid(ip) * wage_occ[io]);
			Enetwage += model.occdist(io) * model.proddist(ip) * netwagegrid.back();
		}
	}

	taxrev = p.corptax * profit - p.lumptransfer;
	for (int io=0; io<nocc; ++io)
		taxrev += p.labtax * wage_occ[io] * labor_occ[io];

	if ( p.taxHHProfitIncome )
		taxrev += p.labtax * p.profdistfracW * profit * (1.0 - p.corptax);
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

	investment = p.depreciation * capital;

	compute_profits();
	compute_factor_prices();
	compute_dividends(p);
	compute_govt(p, model);

	pi = p.pi;
	rnom = rb - p.pi;
	illprice = 1;
	illpricedot = 0;
	illshares = capital + equity_A;
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

	investment = p.depreciation * capital;

	compute_profits();
	compute_factor_prices();
	compute_dividends(p);
	compute_govt(p, model);

	pi = p.pi;
	rnom = rb - p.pi;
	illprice = capital + equity_A;
	illpricedot = 0;
	illshares = initial_equm.illshares;
}

void solve_trans_equilibrium(std::vector<EquilibriumTrans>& trans_equms,
	const Parameters& p, const Model& model, const EquilibriumInitial& initial_equm,
	const EquilibriumFinal& final_equm, const hank_float_type* deltatransvec)
{
	int T = p.Ttransition;
	std::vector<hank_float_type> linv(T);

	for (int it=0; it<T; ++it) {
		double deltatrans = deltatransvec[it];

		trans_equms[it].compute_factors(model);

		trans_equms[it].varieties = trans_equms[it].tfp_N * pow(
			cobb_douglas(trans_equms[it].capital_N, trans_equms[it].labor_N, trans_equms[it].alpha_N),
			trans_equms[it].drs_N);

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
			trans_equms[it].lIK = 0;
		}
		else if ( p.capadjcost > 0 ) {
			trans_equms[it].qinvestment = 0;
			trans_equms[it].invadjust = 0;

			trans_equms[it].lIK = model.capadjcost1inv(trans_equms[it].qcapital - 1.0);

			if (it < T - 1)
				trans_equms[it].qdot = (trans_equms[it+1].qcapital - trans_equms[it].qcapital) / deltatrans;
			else
				trans_equms[it].qdot = 0;

			trans_equms[it].capadjust = model.capadjcost(trans_equms[it].lIK);
			trans_equms[it].ra = (trans_equms[it].rcapital + trans_equms[it].lIK * model.capadjcost1(trans_equms[it].lIK) - trans_equms[it].capadjust + trans_equms[it].qdot)
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

		trans_equms[T-1].lIK = model.capadjcost1inv(trans_equms[T-1].qcapital - 1.0);
		for (int it=0; it<T; ++it)
			trans_equms[it].lIK = linv[it] / trans_equms[it].capital;
	}

	update_equity_variables(p, linv, deltatransvec, final_equm, trans_equms, T);
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

			trans_equms[it].compute_profits();
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
			trans_equms[it].bond = trans_equms[it].equity_B;
		}

		trans_equms[T-1].illpricedot = 0;
		for (int it=T-2; it>=0; --it)
			trans_equms[it].illpricedot = (trans_equms[it+1].illprice - trans_equms[it].illprice) / deltatransvec[it];
	}
}