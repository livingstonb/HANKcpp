#ifndef _ADJUSTMENT_COSTS_H
#define _ADJUSTMENT_COSTS_H

#include <hank_config.h>
#include <functional>
#include <array>
#include <hank.h>

class AdjustmentCosts {
	public:
		AdjustmentCosts() = default;

		AdjustmentCosts(AdjustCostFnRatioMode mode_, bool exponential_costs_,
			double kappa_w_fc_, double kappa_d_fc_, const std::array<double, 5>& kappa_w_,
			const std::array<double, 5>& kappa_d_, double adjcost1max, double dmax);

		std::function<double(double, double)> cost, cost1, cost1inv;

		std::function<double(double)> scale_factor;

		AdjustCostFnRatioMode mode;

		bool exponential_costs;

		double kappa_w_fc, kappa_d_fc, adjcost1max, dmax;

		std::array<double, 5> kappa_w, kappa_d;
};

#endif