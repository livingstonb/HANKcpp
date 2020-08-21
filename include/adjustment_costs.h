#ifndef _ADJUSTMENT_COSTS_H
#define _ADJUSTMENT_COSTS_H

#include <math.h>
#include <hank.h>

class AdjustmentCosts {
	private:
		double cost_fn_exponential(double d, double a) const;

		double cost_fn_other(double d, double a) const;

		double cost_deriv_exponential(double d, double a) const;

		double cost_deriv_other(double d, double a) const;

		double cost_deriv_inv_exponential(double chi, double a) const;

		double cost_deriv_inv_other(double chi, double a) const;

	public:
		AdjustmentCosts() = default;

		AdjustmentCosts(AdjustCostFnRatioMode mode_, bool exponential_costs_,
			double kappa_w_fc_, double kappa_d_fc_, const double kappa_w_[], const double kappa_d_[]);
		
		std::function<double(double, double)> cost;
		std::function<double(double, double)> cost1;
		std::function<double(double, double)> cost1inv;

		double scale_factor(double a) const;

		AdjustCostFnRatioMode mode;
		bool exponential_costs;
		double kappa_w_fc, kappa_d_fc;
		std::array<double, 5> kappa_w;
		std::array<double, 5> kappa_d;
};

#endif