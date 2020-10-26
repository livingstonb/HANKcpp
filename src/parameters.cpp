#include <parameters.h>
#include <algorithm>
#include <vector>
#include <math.h>
#include <string>
#include <utilities.h>

namespace {
	double quadratic_formula(double a, double b, double c);

	double compute_ss_capital_output_ratio(const Parameters& p, double price_W);
}

void Parameters::setup(const Options& opts)
{

	if ( opts.fast ) {
		nb_pos = 10;
		nb_neg = 5;
		na = 10;
	}

	if ( borrowing )
		nb = nb_pos + nb_neg;
	else {
		nb = nb_pos;
		nb_neg = 0;
	}

	nab = na * nb;

	kappa_w[1] = pow(0.5 * (1.0 + kappa_w[2]), -1.0/kappa_w[2]);
	kappa_d[3] = kappa_w[3];
	switch ( depositCostMode ) {
		case DepositCostMode::symmetric:
			// Set kappa_d equal to kappa_w
			kappa_d_fc = kappa_w_fc;
			kappa_d = kappa_w;
			break;
		case DepositCostMode::no_deposit_cost:
			kappa_d_fc = 0;
			kappa_d[0] = 0;
			kappa_d[1] = 100; // This is the optimal depost rate, i.e. higher for lower cost
			kappa_d[2] = 1;
			kappa_d[4] = 0;
			break;
		case DepositCostMode::custom:
			// Leave deposit adj cost parameters at current values
			break;
	}

	if ( illiqWealthTarget.is_mean() )
		targetMeanIllGuess = illiqWealthTarget.value;

	targetMeanIllGuess /= USGDPperHH;
	illiqWealthTarget.value /= USGDPperHH;
	liqWealthTarget.value /= USGDPperHH;

	chi = meanlabeff / (pow(0.7, -riskaver) * pow(hourtarget, 1.0 / frisch));

	update();
}

void Parameters::update() {
	rborr = rb + borrwedge;
	if ( make_profit_correction & (alpha_N == alpha_Y) )
	{
		profdistfracA = alpha_N;
		profdistfracB = 0;
		profdistfracW = 1.0 - alpha_N;
		profdistfracL = 0;
	}

	double price_W = 1.0 - 1.0 / elast;
	target_KY_ratio = compute_ss_capital_output_ratio(*this, price_W);

	if ( global_hank_options->print_diagnostics )
		HANK::print(*this);
}

std::map<std::string, hank_float_type> Parameters::variables_map() const
{
	std::map<std::string, hank_float_type> variables;

	variables.insert({"elast", elast});
	variables.insert({"frisch", frisch});
	variables.insert({"riskaver", riskaver});
	variables.insert({"rho", rho});
	variables.insert({"drs_Y", drs_Y});
	variables.insert({"drs_N", drs_N});
	variables.insert({"alpha_Y", alpha_Y});
	variables.insert({"alpha_N", alpha_N});
	variables.insert({"profdistfracA", profdistfracA});
	variables.insert({"profdistfracB", profdistfracB});
	variables.insert({"profdistfracL", profdistfracL});
	variables.insert({"profdistfracW", profdistfracW});
	variables.insert({"depreciation", depreciation});
	variables.insert({"kappa_w_fc", kappa_w_fc});
	variables.insert({"kappa_w0", kappa_w[0]});
	variables.insert({"kappa_w1", kappa_w[1]});
	variables.insert({"kappa_w2", kappa_w[2]});
	variables.insert({"kappa_w3", kappa_w[3]});
	variables.insert({"kappa_w4", kappa_w[4]});
	variables.insert({"kappa_d_fc", kappa_d_fc});
	variables.insert({"kappa_d0", kappa_d[0]});
	variables.insert({"kappa_d1", kappa_d[1]});
	variables.insert({"kappa_d2", kappa_d[2]});
	variables.insert({"kappa_d4", kappa_d[4]});
	variables.insert({"hourtarget", hourtarget});
	variables.insert({"deathrate", deathrate});
	variables.insert({"rb", rb});
	variables.insert({"rborr", rborr});
	variables.insert({"chi", chi});
	variables.insert({"target_KY_ratio", target_KY_ratio});
	variables.insert({"targetMeanIllGuess", targetMeanIllGuess});
	variables.insert({"drs_Y", drs_Y});
	variables.insert({"drs_Y", drs_Y});
	variables.insert({"drs_Y", drs_Y});
	variables.insert({"drs_Y", drs_Y});

	return variables;
}

namespace HANK {
	void print(const Parameters& p)
	{
		std::map<std::string, hank_float_type> variables = p.variables_map();
		print(variables, "PARAMETERS");
	}
}

namespace {

	double quadratic_formula(double a, double b, double c)
	{
		return (-b + sqrt(pow(b, 2.0) - 4.0 * a * c)) / (2.0 * a);
	}

	double compute_ss_capital_output_ratio(const Parameters& p, double price_W)
	{
		double la, lb, lc;

		la = -p.depreciation;
		lb = p.targetMeanIllGuess * p.depreciation + price_W * p.alpha_Y * p.drs_Y
			+ (1.0 - price_W) * p.alpha_N * p.drs_N
			+ (price_W * (1.0 - p.drs_Y) + (1.0 - price_W) * (1.0 - p.drs_N))
				* (1.0 - p.corptax) * p.profdistfracA;
		lc = -p.targetMeanIllGuess
			* (price_W * p. alpha_Y * p. drs_Y + (1.0 - price_W) * p.alpha_N * p.drs_N);

		return quadratic_formula(la, lb, lc);
	}
}