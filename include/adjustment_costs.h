#ifndef _ADJUSTMENT_COSTS_H
#define _ADJUSTMENT_COSTS_H

#include <math.h>
#include <hank.h>

class AdjustmentCosts {
	public:
		AdjustmentCosts(AdjustCostFnRatioMode mode_, bool exponential_costs_,
			double kappa_w_fc_, double kappa_d_fc_, double kappa_w_[], double kappa_d_[])
			: mode(mode_), exponential_costs(exponential_costs_),
				kappa_w_fc(kappa_w_fc_), kappa_d_fc(kappa_d_fc_) {

			for (int i=0; i<5; ++i) {
				kappa_w[i] = kappa_w_[i];
				kappa_d[i] = kappa_d_[i];
			}
		}

		double cost(double d, double a) {
			double x, scale, acost;

			if ( d == 0 )
				return 0;

			scale = scale_factor(a);
			x = scale * d;

			if ( exponential_costs ) {
				if ( x > 0 )
					acost = kappa_d_fc - x
						+ (1.0 / kappa_d[4]) * (exp(kappa_d[4] * x) - 1.0);
				else
					acost = kappa_w_fc - x
						+ (1.0 / kappa_w[4]) * (exp(kappa_w[4] * x) - 1.0);
			}
			else {
				if ( x > 0 ) {
					acost = kappa_d[0] * fabs(x)
						+ pow(fabs(x), 1.0 + kappa_d[2])
							* pow(kappa_d[1], -kappa_d[2]) / (1.0 + kappa_d[2]);
					acost = kappa_d_fc + scale * acost;
				}
				else {
					acost = kappa_w[0] * fabs(x)
						+ pow(fabs(x), 1.0 + kappa_w[2])
							* pow(kappa_w[1], -kappa_w[2]) / (1.0 + kappa_w[2]);
					acost = kappa_w_fc + scale * acost;
				}
			}

			return acost;
		}

		double cost1(double d, double a) {
			double x, scale, dcost;

			scale = scale_factor(a);
			x = scale * d;

			if ( exponential_costs ) {
				if ( x > 0 )
					dcost = exp(kappa_d[4] * x) - 1.0;
				else
					dcost = exp(kappa_w[4] * x) - 1.0;
			}
			else {
				if ( x > 0 )
					dcost = kappa_d[0] + pow(x / kappa_d[1], kappa_d[2]);
				else
					dcost = -kappa_w[0] - pow(-x / kappa_w[1], kappa_w[2]);
			}
			return dcost;
		}

		double cost1inv(double chi, double a) {return 0.0;}

		double scale_factor(double a) {
			double scale;
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

		AdjustCostFnRatioMode mode;
		bool exponential_costs;
		double kappa_w_fc, kappa_d_fc;
		double kappa_w[5];
		double kappa_d[5];
};

#endif