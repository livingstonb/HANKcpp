#include <steady_state.h>

InitialSteadyState find__initial_steady_state(const Model& model, const Options& options) {
	InitialSteadyState iss = initialize_steady_state(model);

	if ( options.calibrateDiscountRate ) {
		double lrhoL = invlogistic(exp(-0.02));
		double lrhoU = invlogistic(exp(-0.01));
		rtflsp(FnDiscountRate,lrhoL,lrhoU,1.0e-8_8,tolrho,iflag,maxiterrho)
	}
}

InitialSteadyState initialize_steady_state(const Model& model) {
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

	// Compute steady-state target capital-output ratio, given target mean illiq wealth
	double K_totoutput_ratio = compute_ss_capital_output_ratio(price_W,
		p.targetMeanIll, p.depreciation, p.alpha_Y, p.drs_Y, p.alpha_N, p.drs_N);

	// If solving for equilibrium, these are guesses
	double labor_Y = p.hourtarget * p.meanlabeff * model.occdist.dot((1.0-model.occgrid.array()).matrix()) / varieties;
	double labor_N = p.hourtarget * p.meanlabeff * model.occdist.dot(model.occgrid);

	double capital = totoutput * K_totoutput_ratio;
	double rcapital = (price_W * p.alpha_Y * p.drs_Y
		+ (1.0-price_W) * p.alpha_N * p.drs_N) / K_totoutput_ratio;

	// Wholesalers
	double capital_Y = price_W * p.alpha_Y * p.drs_Y * output / rcapital;
	double tfp_Y = output / pow((pow(capital_Y, p.alpha_Y) * pow(labor_Y, 1.0-p.alpha_Y)), p.drs_Y);

	double wage_Y = 0.0;
	if ( (p.alpha_Y  ~= 1.0) & (p.drs_Y ~= 0.0) )
		wage_Y = price_W * (1.0-p.alpha_Y) * p.drs_Y * output / labor_Y;

	double mc_Y = (1.0 / tfp_Y)
		* pow((rcapital/p.alpha_Y), p.alpha_Y)
		* pow((wage_Y / (1.0-p.alpha_Y)), 1.0-p.alpha_Y);

	// Retailers
	double capital_N = grossprofit_R * p.alpha_N * p.drs_N * varieties / rcapital;
	double tfp_N = varieties
		/ pow(pow(capital_N, p.alpha_N)
			* pow(labor_N, 1.0-p.alpha_N), p.drs_N);

	double wage_N = 0.0;
	if ( (p.alpha_N ~= 1.0) & (p.drs_N ~= 0.0) )
		wage_N = grossprofit_R * (1.0 - p.alpha_N) * p.drs_N * varieties / labor_N;

	double mc_N = (1.0/tfp_N) * pow(rcapital/p.alpha_N, p.alpha_N)
		* pow(wage_N/(1.0-p.alpha_N), 1.0-p.alpha_N);

	// Guess labor disutility so that at average wages and average consumption, hours = 1/3
	// (sets C/Y = 0.70) for hand-to-mouth. True hours will be lower.
	double meanwage_y_component = 
	double lmeanwage = wage_Y * model.occdist.dot((1.0-model.occgrid.array()).matrix())
		+ wage_N * model.occdist.dot(model.occgrid);

	double chi = lmeanwage * p.meanlabeff / (pow(0.7, -p.riskaver) * pow(p.hourtarget, 1.0/p.frisch));

	// Output object
	InitialSteadyState ss;

	ss.ra = rcapital - p.depreciation;
	ss.dividend_A = p.profdistfracA * profit * (1.0 - p.corptax);
	ss.equity_A = ss.dividend_A / ss.ra;
	ss.dividend_B = p.profdistfracB * profit * (1.0 - p.corptax);
	ss.equity_B = ss.dividend_B / p.rb;
	ss.wage_Y = wage_Y;
	ss.wage_N = wage_N;
	ss.profit = profit;
	ss.chi = chi;

	ss.netwagegrid = (1.0-p.labtax) * model.yprodgrid.array()
		* (wage_N * model.yoccgrid.array() + wage_Y * (1.0-model.yoccgrid.array()));

	return ss;
}

double compute_ss_capital_output_ratio(double price_W,
	double targetMeanIll, double depreciation, double alpha_Y,
	double drs_Y, double alpha_N, double drs_N)
{
	double la, lb, lc, lKNY;

	la = -depreciation;
	lb = targetMeanIll * depreciation + price_W * alpha_Y * drs_Y
		+ (1.0 - price_W) * alpha_N * drs_N
		+ (price_W * (1.0 - drs_Y) + (1.0 - price_W) * (1.0 - drs_N))
			* profdistfracA;
	lc = -targetMeanIll
		* (price_W * alpha_Y * p.drs_Y + (1.0 - price_W) * alpha_N * drs_N);

	// Quadratic formula
	lKNY = (-lb + sqrt(pow(lb, 2) - 4 * la * lc)) / (2 * la);
	return lKNY;
}