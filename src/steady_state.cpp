#include <steady_state.h>
#include <parameters.h>
#include <model.h>

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
				* p.profdistfracA;
		lc = -p.targetMeanIll
			* (price_W * p. alpha_Y * p. drs_Y + (1.0 - price_W) * p.alpha_N * p.drs_N);

		return quadratic_formula(la, lb, lc);
	}
}

SteadyState::SteadyState(const Model& model_) : model(model_), p(model_.p) {

	totoutput = output * varieties;
	price_W = (1.0 - 1.0 / p.elast);
	grossprofit_W = price_W * output * (1.0 - p.drs_Y);
	netprofit_W = varieties * grossprofit_W;
	grossprofit_R = (1.0 - price_W) * output;
	netprofit_R = varieties * grossprofit_R * (1.0 - p.drs_N);
	profit = netprofit_R + netprofit_W;

	// Compute steady-state target capital-output ratio, given target mean illiq wealth
	K_totoutput_ratio = compute_ss_capital_output_ratio(p, price_W);

	// If solving for equilibrium, these are guesses
	labor_Y = p.hourtarget * p.meanlabeff * model.occdist.dot((1.0-model.occgrid.array()).matrix()) / varieties;
	labor_N = p.hourtarget * p.meanlabeff * model.occdist.dot(model.occgrid);

	// Update other prices
	update();

	// Guess labor disutility so that at average wages and average consumption, hours = 1/3
	// (sets C/Y = 0.70) for hand-to-mouth. True hours will be lower.
	double lmeanwage = wage_Y * model.occdist.dot((1.0-model.occgrid.array()).matrix())
		+ wage_N * model.occdist.dot(model.occgrid);

	chi = lmeanwage * p.meanlabeff / (pow(0.7, -p.riskaver) * pow(p.hourtarget, 1.0/p.frisch));

	yprodgrid = model_.yprodgrid;
}

void SteadyState::update() {
	capital = totoutput * K_totoutput_ratio;
	rcapital = (price_W * p.alpha_Y * p.drs_Y
		+ (1.0-price_W) * p.alpha_N * p.drs_N) / K_totoutput_ratio;

	// Wholesalers
	capital_Y = price_W * p.alpha_Y * p.drs_Y * output / rcapital;
	tfp_Y = output / pow((pow(capital_Y, p.alpha_Y) * pow(labor_Y, 1.0-p.alpha_Y)), p.drs_Y);

	if ( (p.alpha_Y  == 1.0) | (p.drs_Y == 0.0) )
		wage_Y = 0;
	else
		wage_Y = price_W * (1.0 - p.alpha_Y) * p.drs_Y * output / labor_Y;

	mc_Y = marginal_cost(tfp_Y, rcapital, p.alpha_Y, wage_Y);

	// Retailers
	capital_N = grossprofit_R * p.alpha_N * p.drs_N * varieties / rcapital;
	tfp_N = varieties
		/ pow(pow(capital_N, p.alpha_N)
			* pow(labor_N, 1.0-p.alpha_N), p.drs_N);

	if ( (p.alpha_N == 1.0) | (p.drs_N == 0.0) )
		wage_N = 0;
	else
		wage_N = grossprofit_R * (1.0 - p.alpha_N) * p.drs_N * varieties / labor_N;

	mc_Y = marginal_cost(tfp_N, rcapital, p.alpha_N, wage_N);

	ra = rcapital - p.depreciation;
	dividend_A = p.profdistfracA * profit * (1.0 - p.corptax);
	equity_A = dividend_A / ra;
	dividend_B = p.profdistfracB * profit * (1.0 - p.corptax);
	equity_B = dividend_B / p.rb;

	netwagegrid = (1.0 - p.labtax) * model.yprodgrid.array()
		* (wage_N * model.yoccgrid.array() + wage_Y * (1.0 - model.yoccgrid.array()));
}