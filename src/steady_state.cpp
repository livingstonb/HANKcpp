#include <steady_state.h>
#include <parameters.h>
#include <model.h>
#include <hank_eigen_dense.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include <utilities.h>

namespace {

	std::vector<double> repeat(const std::vector<double>& arr_in, int r) {
		std::vector<double> arr_out(arr_in.size() * r);
		int n = arr_in.size();

		auto iter = arr_out.begin();
		for (int i=0; i<n; ++i) {
			std::fill(iter, iter + r, arr_in[i]);
			iter += r;
		}

		return arr_out;
	}
}

SteadyState::SteadyState(const Parameters& p_, const Model& model_) : model(model_), p(p_) {
	price_W = 1.0 - 1.0 / p.elast;
}

void SteadyState::compute(SSType mode) {
	ArrayXr elabshareY, elabshareN, elabfracY, elabfracN;

	assert( static_cast<int>(labor_occ.size()) == p.nocc );

	elabshareY = (1.0 - p.alpha_Y) * price_W * p.drs_Y * model.occYsharegrid;
	elabshareN = (1.0 - p.alpha_N) * (1.0 - price_W) * p.drs_N * model.occNsharegrid;
	labshareY = to_vector(elabshareY);
	labshareN = to_vector(elabshareN);

	capshareY = p.alpha_Y * price_W * p.drs_Y;
	capshareN = p.alpha_N * (1.0 - price_W) * p.drs_N;
	capfracY = capshareY / (capshareY + capshareN);
	capfracN = capshareN / (capshareY + capshareN);
	elabfracY = elabshareY / (elabshareY + elabshareN);
	elabfracN = elabshareN / (elabshareY + elabshareN);

	labfracY = to_vector(elabfracY);
	labfracN = to_vector(elabfracN);
	
	// Output and varieties

	capital_Y = capfracY * capital;
	ArrayXr labor_Yocc = elabfracY * as_eigen<ArrayXr>(labor_occ);
	labor_Y = labor_Yocc.pow(model.occYsharegrid.array()).prod();

	capital_N = capfracN * capital;
	ArrayXr labor_Nocc = elabfracN * as_eigen<ArrayXr>(labor_occ);
	labor_N = labor_Nocc.pow(model.occNsharegrid.array()).prod();
	labor = as_eigen<VectorXr>(labor_occ).dot(model.occdist);

	double yterm = pow(capital_Y, p.alpha_Y) * pow(labor_Y, 1.0-p.alpha_Y);
	double nterm = pow(capital_N, p.alpha_N) * pow(labor_N, 1.0-p.alpha_N);
	if ( mode == SSType::initial ) {
		tfp_Y = output / pow(yterm, p.drs_Y);
		tfp_N = varieties / pow(nterm, p.drs_N);
	}
	else if ( mode == SSType::final ) {
		output = tfp_Y * pow(yterm, p.drs_Y);
		varieties = tfp_N * pow(nterm, p.drs_N);
	}

	investment = p.depreciation * capital;

	compute_profits();
	compute_factor_prices();
	compute_dividends();
	compute_govt();

	if ( mode == SSType::initial ) {
		illprice = 1.0;
		illshares = capital + equity_A;
	}
	else if ( mode == SSType::final ) {
	}
	illpricedot = 0.0;

	if ( global_hank_options->print_diagnostics )
		print_variables();
}

void SteadyState::compute_profits() {
	netprofit_W = price_W * output * (1.0 - p.drs_Y);
	grossprofit_R = (1.0 - price_W) * output / varieties;
	netprofit_R = varieties * (1.0 - p.drs_N) * grossprofit_R;
	profit = netprofit_R + netprofit_W;
}

void SteadyState::compute_factor_prices() {
	double cap_term, lab_term;

	rcapital = (capshareY + capshareN) * output / capital;
	ArrayXr labsharesum = as_eigen<ArrayXr>(labshareN) + as_eigen<ArrayXr>(labshareY);
	wage_occ = to_vector(labsharesum * output * as_eigen<ArrayXr>(labor_occ).inverse());

	// Wholesale
	if ( (p.alpha_Y == 1.0) | (p.drs_Y == 0.0) )
		wage_Y = 0.0;
	else
		wage_Y = price_W * (1.0 - p.alpha_Y) * p.drs_Y * output / labor_Y;

	cap_term = pow(rcapital / p.alpha_Y, p.alpha_Y);
	lab_term = pow(wage_Y / (1.0-p.alpha_Y), 1.0-p.alpha_Y);
	mc_Y = cap_term * lab_term / tfp_Y;

	// Expansion
	if ( (p.alpha_N == 1.0) | (p.drs_N == 0.0) )
		wage_N = 0.0;
	else
		wage_N = grossprofit_R * (1.0 - p.alpha_N) * p.drs_N * varieties / labor_N;

	if ( p.alpha_N > 0.0 ) {
		cap_term = pow(rcapital / p.alpha_N, p.alpha_N);
		lab_term = pow(wage_N / (1.0-p.alpha_N), 1.0-p.alpha_N);
		mc_N = cap_term * lab_term / tfp_N;
	}
	else
		mc_N = wage_N / tfp_N;
}

void SteadyState::compute_dividends() {
	ra = rcapital - p.depreciation;
	dividend_A = p.profdistfracA * profit * (1.0 - p.corptax);
	dividend_B = p.profdistfracB * profit * (1.0 - p.corptax);
	equity_A = dividend_A / ra;
	equity_B = dividend_B / p.rb;
}

void SteadyState::compute_govt() {
	std::vector<double> wage_occ_rep = repeat(wage_occ, model.nprod);
	ArrayXr wage_occ_ygrid = as_eigen<ArrayXr>(wage_occ_rep);
	VectorXr enetwagegrid = (1.0 - p.labtax) * model.yprodgrid.array() * wage_occ_ygrid;
	netwagegrid = to_vector(enetwagegrid);
	Enetwage = enetwagegrid.dot(model.ydist);

	taxrev = p.labtax * as_eigen<VectorXr>(wage_occ).dot(as_eigen<VectorXr>(labor_occ))
		- p.lumptransfer + p.corptax * profit;

	if ( p.taxHHProfitIncome )
		taxrev += p.labtax * p.profdistfracW * profit * (1.0 - p.corptax);
}

void SteadyState::print_variables() const {
	std::cout << '\n';
	horzline();
	std::cout << "COMPUTED VALUES, STEADY STATE:\n";

	std::vector<std::string> names;
	std::vector<double> values;

	names.push_back("price_W");
	values.push_back(price_W);

	names.push_back("capshareY");
	values.push_back(capshareY);

	names.push_back("capshareN");
	values.push_back(capshareN);

	names.push_back("capital_Y");
	values.push_back(capital_Y);

	names.push_back("capital_N");
	values.push_back(capital_N);

	names.push_back("labor_Y");
	values.push_back(labor_Y);

	names.push_back("labor_N");
	values.push_back(labor_N);

	names.push_back("labor");
	values.push_back(labor);

	names.push_back("profit");
	values.push_back(profit);

	names.push_back("ra");
	values.push_back(ra);

	names.push_back("dividend_A");
	values.push_back(dividend_A);

	names.push_back("dividend_B");
	values.push_back(dividend_B);

	names.push_back("equity_A");
	values.push_back(equity_A);

	names.push_back("equity_B");
	values.push_back(equity_B);

	names.push_back("E[netwage]");
	values.push_back(Enetwage);

	print_values(names, values);

	horzline();
}