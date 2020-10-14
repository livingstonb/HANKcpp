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

void Parameters::setup(const Options& opts) {

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
	if ( make_profit_correction & (alpha_N == alpha_Y) ) {
		profdistfracA = alpha_N;
		profdistfracB = 0;
		profdistfracW = 1.0 - alpha_N;
		profdistfracL = 0;
	}

	double price_W = 1.0 - 1.0 / elast;
	target_KY_ratio = compute_ss_capital_output_ratio(*this, price_W);

	if ( global_hank_options->print_diagnostics )
		print_variables();
}

void Parameters::print_variables() const {
	std::cout << '\n';
	HankUtilities::horzline();
	std::cout << "PARAMETER VALUES:\n";

	std::vector<std::string> names;
	std::vector<double> values;

	names.push_back("elast");
	values.push_back(elast);

	names.push_back("frisch");
	values.push_back(frisch);

	names.push_back("riskaver");
	values.push_back(riskaver);

	names.push_back("rho");
	values.push_back(rho);

	names.push_back("drs_Y");
	values.push_back(drs_Y);

	names.push_back("drs_N");
	values.push_back(drs_N);

	names.push_back("alpha_Y");
	values.push_back(alpha_Y);

	names.push_back("alpha_N");
	values.push_back(alpha_N);

	names.push_back("depreciation");
	values.push_back(depreciation);

	names.push_back("kappa_w_fc");
	values.push_back(kappa_w_fc);

	for (int i=0; i<5; ++i) {
		names.push_back("kappa_w"+std::to_string(i));
		values.push_back(kappa_w[i]);
	}

	names.push_back("kappa_d_fc");
	values.push_back(kappa_d_fc);

	for (int i=0; i<5; ++i) {
		names.push_back("kappa_d"+std::to_string(i));
		values.push_back(kappa_d[i]);
	}

	names.push_back("hourtarget");
	values.push_back(hourtarget);

	names.push_back("deathrate");
	values.push_back(deathrate);

	names.push_back("rb");
	values.push_back(rb);

	names.push_back("rborr");
	values.push_back(rborr);

	names.push_back("chi");
	values.push_back(chi);

	names.push_back("target_KY_ratio");
	values.push_back(target_KY_ratio);

	names.push_back("targetMeanIllGuess");
	values.push_back(targetMeanIllGuess);

	HankUtilities::print_values(names, values);

	HankUtilities::horzline();
}

namespace {

	double quadratic_formula(double a, double b, double c) {
		return (-b + sqrt(pow(b, 2.0) - 4.0 * a * c)) / (2.0 * a);
	}

	double compute_ss_capital_output_ratio(const Parameters& p, double price_W) {
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