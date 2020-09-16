#include <steady_state.h>
#include <parameters.h>
#include <model.h>
#include <hank_eigen_dense.h>
#include <math.h>
#include <algorithm>
#include <iostream>

namespace {
	double marginal_cost(double tfp, double r, double alpha, double wage) {
		return (1.0 / tfp) * pow(r / alpha, alpha) * pow(wage / (1.0 - alpha), 1.0 - alpha);
	}

	double quadratic_formula(double a, double b, double c) {
		return (-b + sqrt(pow(b, 2.0) - 4.0 * a * c)) / (2.0 * a);
	}

	double compute_ss_capital_output_ratio(const Parameters& p, double price_W) {
		double la, lb, lc;

		la = -p.depreciation;
		lb = p.targetMeanIll * p.depreciation + price_W * p.alpha_Y * p.drs_Y
			+ (1.0 - price_W) * p.alpha_N * p.drs_N
			+ (price_W * (1.0 - p.drs_Y) + (1.0 - price_W) * (1.0 - p.drs_N))
				* (1.0 - p.corptax) * p.profdistfracA;
		lc = -p.targetMeanIll
			* (price_W * p. alpha_Y * p. drs_Y + (1.0 - price_W) * p.alpha_N * p.drs_N);

		return quadratic_formula(la, lb, lc);
	}

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
	target_KY_ratio = compute_ss_capital_output_ratio(p, price_W);
}

void SteadyState::set(const double x[], SSType mode) {
	if ( mode == SSType::initial ) {
		// Guess rho and labor occ, and chi
		labor_occ.clear();
		for (int io=0; io<p.nocc; ++io)
			labor_occ.push_back(x[io+1]);

		capital = x[p.nocc+1];
	}
	else if ( mode == SSType::final ) {
	}
}

void SteadyState::compute(SSType mode) {
	ArrayXd elabshareY, elabshareN, elabfracY, elabfracN;

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
	ArrayXd labor_Yocc = elabfracY * as_eigen<ArrayXd>(labor_occ);
	labor_Y = labor_Yocc.pow(model.occYsharegrid.array()).prod();

	capital_N = capfracN * capital;
	ArrayXd labor_Nocc = elabfracN * as_eigen<ArrayXd>(labor_occ);
	labor_N = labor_Nocc.pow(model.occNsharegrid.array()).prod();
	labor = as_eigen<VectorXd>(labor_occ).dot(model.occdist);

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
	double_array_vector labsharesum = as_eigen<ArrayXd>(labshareN) + as_eigen<ArrayXd>(labshareY);
	wage_occ = to_vector(labsharesum * output * as_eigen<ArrayXd>(labor_occ).inverse());

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
	ArrayXd wage_occ_ygrid = as_eigen<ArrayXd>(wage_occ_rep);
	double_vector enetwagegrid = (1.0 - p.labtax) * model.yprodgrid.array() * wage_occ_ygrid;
	netwagegrid = to_vector(enetwagegrid);

	taxrev = p.labtax * as_eigen<VectorXd>(wage_occ).dot(as_eigen<VectorXd>(labor_occ))
		- p.lumptransfer + p.corptax * profit;

	if ( p.taxHHProfitIncome )
		taxrev += p.labtax * p.profdistfracW * profit * (1.0 - p.corptax);
}

