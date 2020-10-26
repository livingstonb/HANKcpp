#include <adjustment_costs.h>
#include <math.h>
#include <iostream>
#include <functional>

using namespace std::placeholders;

namespace {
	double cost_fn_exponential(AdjustmentCosts* adjcosts, double d, double a);

	double cost_fn_other(AdjustmentCosts* adjcosts, double d, double a);

	double cost_deriv_exponential(AdjustmentCosts* adjcosts, double d, double a);

	double cost_deriv_other(AdjustmentCosts* adjcosts, double d, double a);

	double cost_deriv_inv_exponential(AdjustmentCosts* adjcosts, double chi, double a);

	double cost_deriv_inv_other(AdjustmentCosts* adjcosts, double chi, double a);
}

AdjustmentCosts::AdjustmentCosts(AdjustCostFnRatioMode mode_, bool exponential_costs_,
	double kappa_w_fc_, double kappa_d_fc_, const std::array<double, 5>& kappa_w_,
	const std::array<double, 5>& kappa_d_, double adjcost1max_, double dmax_)
	: mode(mode_), exponential_costs(exponential_costs_),
		kappa_w_fc(kappa_w_fc_), kappa_d_fc(kappa_d_fc_), kappa_w(kappa_w_), kappa_d(kappa_d_) {

	dmax = dmax_;
	adjcost1max = adjcost1max_;

	double kappaw3 = kappa_w[3];
	switch ( mode ) {
		case AdjustCostFnRatioMode::none:
			scale_factor = [](double) {return 1.0;};
			break;
		case AdjustCostFnRatioMode::linear:
			scale_factor = [=](double a) {return 1.0 / (kappaw3 + a);};
			break;
		case AdjustCostFnRatioMode::max:
			scale_factor = [=](double a) {return 1.0 / fmax(a, kappaw3);};
			break;
	}

	if ( exponential_costs ) {
		cost = std::bind(cost_fn_exponential, this, _1, _2);
		cost1 = std::bind(cost_deriv_exponential, this, _1, _2);
		cost1inv = std::bind(cost_deriv_inv_exponential, this, _1, _2);
	}
	else {
		cost = std::bind(cost_fn_other, this, _1, _2);
		cost1 = std::bind(cost_deriv_other, this, _1, _2);
		cost1inv = std::bind(cost_deriv_inv_other, this, _1, _2);
	}
}

namespace {
	double cost_fn_exponential(AdjustmentCosts* adjcosts, double d, double a)
	{
		double x = adjcosts->scale_factor(a) * d;

		if ( x == 0 )
			return 0;
		else if ( x > 0 )
			return adjcosts->kappa_d_fc - x + (1.0 / adjcosts->kappa_d[4]) * (exp(adjcosts->kappa_d[4] * x) - 1.0);
		else
			return adjcosts->kappa_w_fc - x + (1.0 / adjcosts->kappa_w[4]) * (exp(adjcosts->kappa_w[4] * x) - 1.0);
	}

	double cost_fn_other(AdjustmentCosts* adjcosts, double d, double a)
	{
		double scale = adjcosts->scale_factor(a);
		double fcost, x = scale * d;
		const double* kappa = nullptr;

		if ( x == 0 )
			return 0;
		else if ( x > 0 ) {
			kappa = adjcosts->kappa_d.data();
			fcost = adjcosts->kappa_d_fc;
		}
		else {
			kappa = adjcosts->kappa_w.data();
			fcost = adjcosts->kappa_w_fc;
		}
		return fcost + (kappa[0] * fabs(x) + pow(fabs(x), 1.0 + kappa[2])
								* pow(kappa[1], -kappa[2]) / (1.0 + kappa[2])) / scale;
	}

	double cost_deriv_exponential(AdjustmentCosts* adjcosts, double d, double a)
	{
		double x = adjcosts->scale_factor(a) * d;

		if ( x == 0 )
			return 0;
		else if ( x > 0 )
			return exp(adjcosts->kappa_d[4] * x) - 1.0;
		else
			return exp(adjcosts->kappa_w[4] * x) - 1.0;
	}

	double cost_deriv_other(AdjustmentCosts* adjcosts, double d, double a)
	{
		double x = adjcosts->scale_factor(a) * d;

		if ( d == 0 )
			return 0;
		else if ( d > 0 )
			return adjcosts->kappa_d[0] + pow(x / adjcosts->kappa_d[1], adjcosts->kappa_d[2]);
		else
			return -adjcosts->kappa_w[0] - pow(-x / adjcosts->kappa_w[1], adjcosts->kappa_w[2]);
	}

	double cost_deriv_inv_exponential(AdjustmentCosts* adjcosts, double chi, double a)
	{
		double scale = adjcosts->scale_factor(a);
		
		if ( chi > adjcosts->adjcost1max )
			return adjcosts->dmax / scale;
		else if ( chi < -adjcosts->adjcost1max )
			return -adjcosts->dmax / scale;

		if ( chi == 0 )
			return 0;
		else if ( chi > 0 )
			return (1 / adjcosts->kappa_d[4]) * log(1 + chi) / scale;
		else
			return (1 / adjcosts->kappa_w[4]) * log(1 + chi) / scale;
	}

	double cost_deriv_inv_other(AdjustmentCosts* adjcosts, double chi, double a)
	{
		double scale = adjcosts->scale_factor(a);

		if ( chi > adjcosts->adjcost1max )
			return adjcosts->dmax / scale;
		else if ( chi < -adjcosts->adjcost1max )
			return -adjcosts->dmax / scale;

		if ( (chi >= -adjcosts->kappa_w[0]) & (chi <= adjcosts->kappa_d[0]) )
			return 0;
		else if ( chi > adjcosts->kappa_d[0] )
			return adjcosts->kappa_d[1] * pow(chi - adjcosts->kappa_d[0], 1 / adjcosts->kappa_d[2]) / scale;
		else
			return -adjcosts->kappa_w[1] * pow(-chi - adjcosts->kappa_w[0], 1 /adjcosts->kappa_w[2]) / scale;
	}
}