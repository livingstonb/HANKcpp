#include <parameters.h>
#include <algorithm>
#include <vector>
#include <math.h>

namespace {

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
}

void Parameters::setup(const Options& opts) {

	if ( opts.fast ) {
		nb_pos = 10;
		nb_neg = 5;
		na = 10;
	}

	if ( borrowing )
		nb = nb_pos + nb_neg;
	else
		nb = nb_pos;

	ny = nprod * nocc;
	naby = na * nb * ny;
	nab = na * nb;

	rborr = rb + borrwedge;

	kappa_w[1] = pow(0.5 * (1.0 + kappa_w[2]), -1.0/kappa_w[2]);
	kappa_d[3] = kappa_w[3];
	switch ( opts.depositCostMode ) {
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

	double price_W = 1.0 - 1.0 / elast;
	target_KY_ratio = compute_ss_capital_output_ratio(*this, price_W);
	chi = meanlabeff / (pow(0.7, -riskaver) * pow(hourtarget, 1.0 / frisch));
}

void Parameters::update() {
	rborr = rb + borrwedge;
	target_KY_ratio = compute_ss_capital_output_ratio(*this, 1.0 - 1.0 / elast);
}