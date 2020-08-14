#include <steady_state.h>

void solve_initial_steady_state(const Model& model, const Options& options)
{
	std::cout << "Solving for initial steady state" << '\n';

	if ( (~options.equilibriumR) & (~options.calibrateDiscountRate) ) {
			solve_initial_prices(model);
			// iterate_bellman(model);
		}
}

void solve_initial_prices(const Model& model) {
	bool initialSS = true;
	double lmeanwage;
	const Parameters& p = model.p;

	// Normalize steady state Y = N = P = P_R = 1
	double output = 1.0;
	double varieties = 1.0;
	double totoutput = output * varieties;

	// Steady-state prices and profits determined by elasticity of substitution
	double price_W = (1.0 - 1.0 / p.elast);
	double grossprofit_W = price_W * output * (1.0 - p.drs_Y);
	double netprofit_W = varieties * grossprofit_W;
	double grossprofit_R = (1.0 - price_W) * output;
	double netprofit_R = varieties * grossprofit_R * (1.0 - p.drs_N);
	double profit = netprofit_R + netprofit_W;

	// Compute steady-state target capital-output ratio
	double K_totoutput_ratio = compute_ss_capital_output_ratio(p, price_W);

	// If solving for equilibrium, these are guesses
	double labor_Y = p.hourtarget * p.meanlabeff * model.occdist.dot((1.0-model.occgrid.array()).matrix()) / varieties;
	double labor_N = p.hourtarget * p.meanlabeff * model.occdist.dot(model.occgrid);

	double capital = totoutput * K_totoutput_ratio;
	double rcapital = (price_W * p.alpha_Y * p.drs_Y
		+ (1.0-price_W) * p.alpha_N * p.drs_N) / K_totoutput_ratio;

	double capital_Y = price_W * p.alpha_Y * p.drs_Y * output / rcapital;
	double tfp_Y = output / pow((pow(capital_Y, p.alpha_Y) * pow(labor_Y, 1.0-p.alpha_Y)), p.drs_Y);

	double wage_Y = price_W * (1.0-p.alpha_Y) * p.drs_Y * output / labor_Y;
	if ( (p.alpha_Y == 1.0) | (p.drs_Y == 0.0) )
		wage_Y = 0.0;

	double mc_Y = (1.0 / tfp_Y)
		* pow((rcapital/p.alpha_Y), p.alpha_Y)
		* pow((wage_Y / (1.0-p.alpha_Y)), 1.0-p.alpha_Y);

	double capital_N = grossprofit_R * p.alpha_N * p.drs_N * varieties / rcapital;
	double tfp_N = varieties
		/ pow(pow(capital_N, p.alpha_N)
			* pow(labor_N, 1.0-p.alpha_N), p.drs_N);

	double wage_N = grossprofit_R * (1.0 - p.alpha_N) * p.drs_N * varieties / labor_N;
	if ( (p.alpha_N == 1.0) | (p.drs_N == 0.0) )
		wage_N = 0.0;

	double mc_N = (1.0/tfp_N) * pow(rcapital/p.alpha_N, p.alpha_N)
		* pow(wage_N/(1.0-p.alpha_N), 1.0-p.alpha_N);

	double ra = rcapital - p.depreciation;
	double dividend_A = p.profdistfracA * profit * (1.0 - p.corptax);
	double equity_A = dividend_A / ra;
	double dividend_B = p.profdistfracB * profit * (1.0 - p.corptax);
	double equity_B = dividend_B / p.rb;

	double_vector netwagegrid = (1.0-p.labtax) * model.yprodgrid.array()
		* (wage_N * model.yoccgrid.array() + wage_Y * (1.0-model.yoccgrid.array()));

}

double compute_ss_capital_output_ratio(
	const Parameters& p, double lprice_W)
{
	double la, lb, lc, lKNY;

	la = -p.depreciation;
	lb = p.targetMeanIll * p.depreciation + lprice_W * p.alpha_Y * p.drs_Y
		+ (1.0 - lprice_W) * p.alpha_N * p.drs_N
		+ (lprice_W * (1.0 - p.drs_Y) + (1.0 - lprice_W) * (1.0 - p.drs_N))
			* p.profdistfracA;
	lc = -p.targetMeanIll
		* (lprice_W * p.alpha_Y * p.drs_Y + (1.0 - lprice_W) * p.alpha_N * p.drs_N);
	lKNY = (-lb + sqrt(pow(lb, 2) - 4 * la * lc)) / (2 * la);

	return lKNY;
}