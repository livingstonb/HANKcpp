#include <initial_steady_state.h>

void solve_initial_steady_state(const Model& model)
{
	std::cout << "Solving for initial steady state" << '\n';
	solve_initial_prices(model);
}

void solve_initial_prices(const Model& model) {
	bool initialSS = true;
	double lmeanwage;

	// Normalize steady state Y = N = P= P_R = 1
	double output = 1.0;
	double varieties = 1.0;
	double totoutput = output * varieties;

	// Steady-state prices and profits determined by elasticity of substitution
	double price_W = (1.0 - 1.0 / model.p.elast);
	double grossprofit_W = price_W * output * (1.0 - model.p.drs_Y);
	double netprofit_W = varieties * grossprofit_W;
	double grossprofit_R = (1.0 - price_W) * output;
	double netprofit_R = varieties * grossprofit_R * (1.0 - model.p.drs_N);
	double profit = netprofit_R + netprofit_W;

	double K_totoutput_ratio = compute_ss_capital_output_ratio(model.p, price_W);
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