#include <adjustment_costs.h>
#include <math.h>
#include <hank.h>

AdjustmentCosts::AdjustmentCosts(AdjustCostFnRatioMode mode_, bool exponential_costs_,
	double kappa_w_fc_, double kappa_d_fc_, const double kappa_w_[], const double kappa_d_[])
	: mode(mode_), exponential_costs(exponential_costs_),
		kappa_w_fc(kappa_w_fc_), kappa_d_fc(kappa_d_fc_) {

	for (int i=0; i<5; ++i) {
		kappa_w[i] = kappa_w_[i];
		kappa_d[i] = kappa_d_[i];
	}

	if ( exponential_costs ) {
		cost = [this](double d, double a) {return cost_fn_exponential(d, a);};
		cost1 = [this](double d, double a) {return cost_deriv_exponential(d, a);};
		cost1inv = [this](double chi, double a) {return cost_deriv_inv_exponential(chi, a);};
	}
	else {
		cost = [this](double d, double a) {return cost_fn_other(d, a);};
		cost1 = [this](double d, double a) {return cost_deriv_other(d, a);};
		cost1inv = [this](double chi, double a) {return cost_deriv_inv_other(chi, a);};
	}
}

double AdjustmentCosts::scale_factor(double a) const {
	double scale = 0.0;
	switch (mode) {
		case AdjustCostFnRatioMode::none:
			scale = 1;
			break;
		case AdjustCostFnRatioMode::linear:
			scale = 1 / (kappa_w[3] + a);
			break;
		case AdjustCostFnRatioMode::max:
			scale = 1 / fmax(a, kappa_w[3]);
			break;
	}
	return scale;
}

double AdjustmentCosts::cost_fn_exponential(double d, double a) const {
	double x = scale_factor(a) * d;

	if ( x == 0 )
		return 0;
	else if ( x > 0 )
		return kappa_d_fc - x + (1.0 / kappa_d[4]) * (exp(kappa_d[4] * x) - 1.0);
	else
		return kappa_w_fc - x + (1.0 / kappa_w[4]) * (exp(kappa_w[4] * x) - 1.0);
}

double AdjustmentCosts::cost_fn_other(double d, double a) const {
	double scale = scale_factor(a);
	double fcost, x = scale * d;
	const double* kappa = NULL;

	if ( x == 0 )
		return 0;
	else if ( x > 0 ) {
		kappa = kappa_d.data();
		fcost = kappa_d_fc;
	}
	else {
		kappa = kappa_w.data();
		fcost = kappa_w_fc;
	}
	return fcost + scale * (kappa[0] * fabs(x) + pow(fabs(x), 1.0 + kappa[2])
							* pow(kappa[1], -kappa[2]) / (1.0 + kappa[2]));
}

double AdjustmentCosts::cost_deriv_exponential(double d, double a) const {
	double x = scale_factor(a) * d;

	if ( x == 0 )
		return 0;
	else if ( x > 0 )
		return exp(kappa_d[4] * x) - 1.0;
	else
		return exp(kappa_w[4] * x) - 1.0;
}

double AdjustmentCosts::cost_deriv_other(double d, double a) const {
	double x = scale_factor(a) * d;

	if ( x == 0 )
		return 0;
	else if ( x > 0 )
		return kappa_d[0] + pow(x / kappa_d[1], kappa_d[2]);
	else
		return -kappa_w[0] - pow(-x / kappa_w[1], kappa_w[2]);
}

double AdjustmentCosts::cost_deriv_inv_exponential(double chi, double a) const {
	double scale = scale_factor(a);
	if ( chi == 0 )
		return 0;
	else if ( chi > 0 )
		return (1 / kappa_d[4]) * log(1 + chi) / scale;
	else
		return (1 / kappa_w[4]) * log(1 + chi) / scale;
}

double AdjustmentCosts::cost_deriv_inv_other(double chi, double a) const {
	double scale = scale_factor(a);
	if ( (chi >= -kappa_w[0]) & (chi <= kappa_d[0]) )
		return 0;
	else if ( chi > kappa_d[0] )
		return kappa_d[1] * pow(chi - kappa_d[0], 1 / kappa_d[2]) / scale;
	else
		return -kappa_w[1] * pow(-chi - kappa_w[0], 1 /kappa_w[2]) / scale;
}