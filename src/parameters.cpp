#include <parameters.h>
#include <algorithm>
#include <vector>
#include <math.h>
#include <string>

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

	void print_parameter(const std::string& pname, double value, bool insert_endline) {
		std::cout << "  " << pname << " = " << value;

		if ( insert_endline )
			std::cout << '\n';
	}

	void print_parameter(const std::string& pname, double value) {
		print_parameter(pname, value, true);
	}
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

	rborr = rb + borrwedge;

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

	double price_W = 1.0 - 1.0 / elast;
	target_KY_ratio = compute_ss_capital_output_ratio(*this, price_W);
	chi = meanlabeff / (pow(0.7, -riskaver) * pow(hourtarget, 1.0 / frisch));
}

void Parameters::update() {
	rborr = rb + borrwedge;
	target_KY_ratio = compute_ss_capital_output_ratio(*this, 1.0 - 1.0 / elast);

	if ( global_hank_options->print_diagnostics )
		print_values();
}

void Parameters::print_values() const {
	horzline();

	std::cout << "PARAMETER VALUES:\n";
	print_parameter("elast", elast);
	print_parameter("frisch", frisch);
	print_parameter("riskaver", riskaver);
	print_parameter("rho", rho);
	print_parameter("drs_Y", drs_Y);
	print_parameter("drs_N", drs_N);
	print_parameter("alpha_Y", alpha_Y);
	print_parameter("alpha_N", alpha_N);
	print_parameter("depreciation", depreciation);

	print_parameter("kappa_w_fc", kappa_w_fc);
	for (int i=0; i<5; ++i)
		print_parameter("kappa_w"+std::to_string(i), kappa_w[i]);

	print_parameter("kappa_dw_fc", kappa_d_fc);
	for (int i=0; i<5; ++i)
		print_parameter("kappa_d"+std::to_string(i), kappa_d[i]);

	print_parameter("hourtarget", hourtarget);
	print_parameter("deathrate", deathrate);
	print_parameter("rb", rb);
	print_parameter("rborr", rborr);
	print_parameter("chi", chi);
	print_parameter("target_KY_ratio", target_KY_ratio);
	print_parameter("targetMeanIllGuess", targetMeanIllGuess, false);

	horzline();
	std::cout << '\n';
}