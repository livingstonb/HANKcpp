#include <equilibrium.h>
#include <parameters.h>
#include <model.h>
#include <algorithm>

namespace
{
	double cobb_douglas(hank_float_type capital, hank_float_type labor, hank_float_type alpha);
}

void EquilibriumElement::create_initial_steady_state(const Parameters& p, const Model& model)
{
	set_parameters(p);
	nprod = model.nprod;

	if ( labor_occ.size() == 0 )
		for (int io=0; io<nocc; ++io)
			labor_occ.push_back(p.hourtarget * p.meanlabeff * model.occdist[io]);

	if ( capital == CapitalNotSet )
		capital = p.target_KY_ratio;

	compute_factors(model);

	tfp_Y = output / pow(cobb_douglas(capital_Y, labor_Y, alpha_Y), drs_Y);
	tfp_N = varieties / pow(cobb_douglas(capital_N, labor_N, alpha_N), drs_N);

	investment = p.depreciation * capital;

	compute_profits();
	compute_factor_prices();
	compute_dividends(p);
	compute_govt(p, model);
}

void EquilibriumElement::set_parameters(const Parameters& p)
{
	alpha_Y = p.alpha_Y;
	alpha_N = p.alpha_N;
	price_W = 1.0 - 1.0 / p.elast;
	drs_Y = p.drs_Y;
	drs_N = p.drs_N;
	nocc = p.nocc;
}

void EquilibriumElement::compute_factors(const Model& model)
{
	// Capital
	capshareY = alpha_Y * price_W * drs_Y;
	capshareN = alpha_N * (1.0 - price_W) * drs_N;
	capfracY = capshareY / (capshareY + capshareN);
	capfracN = capshareN / (capshareY + capshareN);
	capital_Y = capfracY * capital;
	capital_N = capfracN * capital;

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
}

void EquilibriumElement::compute_profits()
{
	netprofit_W = price_W * output * (1.0 - drs_Y);
	grossprofit_R = (1.0 - price_W) * output / varieties;
	netprofit_R = varieties * (1.0 - drs_N) * grossprofit_R;
	profit = netprofit_R + netprofit_W;
}

void EquilibriumElement::compute_factor_prices()
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

void EquilibriumElement::compute_dividends(const Parameters& p)
{
	ra = rcapital - p.depreciation;
	rb = p.rb;
	dividend_A = p.profdistfracA * profit * (1.0 - p.corptax);
	dividend_B = p.profdistfracB * profit * (1.0 - p.corptax);
	equity_A = dividend_A / ra;
	equity_B = dividend_B / p.rb;
}

void EquilibriumElement::compute_govt(const Parameters& p, const Model& model)
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

namespace
{
	double cobb_douglas(hank_float_type capital, hank_float_type labor, hank_float_type alpha)
	{
		return pow(capital, alpha) * pow(labor, 1.0 - alpha);
	}
}