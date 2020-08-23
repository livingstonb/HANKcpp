#ifndef _ADJUSTMENT_COSTS_H
#define _ADJUSTMENT_COSTS_H

#include <hank_config.h>
#include <functional>
#include <array>
#include <hank_types.h>

enum class AdjustCostFnRatioMode;

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
			double kappa_w_fc_, double kappa_d_fc_, const std::array<double, 5>& kappa_w_,
			const std::array<double, 5>& kappa_d_);
		
		double scale_factor(double a) const;

		std::function<double(double, double)> cost, cost1, cost1inv;

		AdjustCostFnRatioMode mode;

		bool exponential_costs;

		double kappa_w_fc, kappa_d_fc;

		std::array<double, 5> kappa_w, kappa_d;
};

#endif