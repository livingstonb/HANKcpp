#include <initial_steady_state.h>

void solve_initial_steady_state(const Model& model)
{
	std::cout << "Solving for initial steady state" << '\n';
	solve_initial_prices(model);
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
	double labor_Y = model.hourtarget * model.meanlabeff * model.occdist.dot(1.0-model.occgrid) / varieties;
	double labor_N = model.hourtarget * model.meanlabeff * model.occdist.dot(model.occgrid);

	double capital = totoutput * K_totoutput_ratio;
	double rcapital = (price_W * model.alpha_Y * model.drs_Y
		+ (1.0-price_W) * model.alpha_N * model.drs_N) / K_totoutput_ratio;
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