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

	// Normalize steady state Y = N = P = P_R = 1
	double output = 1.0;
	double varieties = 1.0;
	double totoutput = output * varieties;

	// Steady-state prices and profits determined by elasticity of substitution
	double price_W = (1.0 - 1.0 / model.p.elast);
	double grossprofit_W = price_W * output * (1.0 - model.drs_Y);
	double netprofit_W = varieties * grossprofit_W;
	double grossprofit_R = (1.0 - price_W) * output;
	double netprofit_R = varieties * grossprofit_R * (1.0 - model.drs_N);
	double profit = netprofit_R + netprofit_W;

	// Compute steady-state target capital-output ratio
	double K_totoutput_ratio = compute_ss_capital_output_ratio(model, price_W);

	// If solving for equilibrium, these are guesses
	double labor_Y = model.hourtarget * model.meanlabeff * model.occdist.dot((1.0-model.occgrid.array()).matrix()) / varieties;
	double labor_N = model.hourtarget * model.meanlabeff * model.occdist.dot(model.occgrid);

	double capital = totoutput * K_totoutput_ratio;
	double rcapital = (price_W * model.alpha_Y * model.drs_Y
		+ (1.0-price_W) * model.alpha_N * model.drs_N) / K_totoutput_ratio;

	double capital_Y = price_W * model.alpha_Y * model.drs_Y * output / rcapital;
	double tfp_Y = output / pow((pow(capital_Y, model.alpha_Y) * pow(labor_Y, 1.0-model.alpha_Y)), model.drs_Y);

	double wage_Y = price_W * (1.0-model.alpha_Y) * model.drs_Y * output / labor_Y;
	if ( (model.alpha_Y == 1.0) | (model.drs_Y == 0.0) )
		wage_Y = 0.0;

	double mc_Y = (1.0 / tfp_Y)
		* pow((rcapital/model.alpha_Y), model.alpha_Y)
		* pow((wage_Y / (1.0-model.alpha_Y)), 1.0-model.alpha_Y);

	double capital_N = grossprofit_R * model.alpha_N * model.drs_N * varieties / rcapital;
	double tfp_N = varieties
		/ pow(pow(capital_N, model.alpha_N)
			* pow(labor_N, 1.0-model.alpha_N), model.drs_N);

	double wage_N = grossprofit_R * (1.0 - model.alpha_N) * model.drs_N * varieties / labor_N;
	if ( (model.alpha_N == 1.0) | (model.drs_N == 0.0) )
		wage_N = 0.0;

	double mc_N = (1.0/tfp_N) * pow(rcapital/model.alpha_N, model.alpha_N)
		* pow(wage_N/(1.0-model.alpha_N), 1.0-model.alpha_N);

	double ra = rcapital - model.depreciation;
	double dividend_A = model.profdistfracA * profit * (1.0 - model.corptax);
	double equity_A = dividend_A / ra;
	double dividend_B = model.profdistfracB * profit * (1.0 - model.corptax);
	double equity_B = dividend_B / model.rb;

	double_vector netwagegrid = (1.0-model.labtax) * model.yprodgrid.array()
		* (wage_N * model.yoccgrid.array() + wage_Y * (1.0-model.yoccgrid.array()));

}

double compute_ss_capital_output_ratio(
	const Model& model, double lprice_W)
{
	double la, lb, lc, lKNY;

	la = -model.depreciation;
	lb = model.targetMeanIll * model.depreciation + lprice_W * model.alpha_Y * model.drs_Y
		+ (1.0 - lprice_W) * model.alpha_N * model.drs_N
		+ (lprice_W * (1.0 - model.drs_Y) + (1.0 - lprice_W) * (1.0 - model.drs_N))
			* model.profdistfracA;
	lc = -model.targetMeanIll
		* (lprice_W * model.alpha_Y * model.drs_Y + (1.0 - lprice_W) * model.alpha_N * model.drs_N);
	lKNY = (-lb + sqrt(pow(lb, 2) - 4 * la * lc)) / (2 * la);

	return lKNY;
}