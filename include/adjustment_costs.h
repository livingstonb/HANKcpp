#ifndef _ADJUSTMENT_COSTS_H
#define _ADJUSTMENT_COSTS_H

#include <math.h>
#include <hank.h>

class AdjustmentCosts {
	private:
		double cost_fn_exponential(double d, double a) const {
			double x = scale_factor(a) * d;

			if ( x == 0 )
				return 0;
			else if ( x > 0 )
				return kappa_d_fc - x + (1.0 / kappa_d[4]) * (exp(kappa_d[4] * x) - 1.0);
			else
				return kappa_w_fc - x + (1.0 / kappa_w[4]) * (exp(kappa_w[4] * x) - 1.0);
		}

		double cost_fn_other(double d, double a) const {
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

		double cost_deriv_exponential(double d, double a) const {
			double x = scale_factor(a) * d;

			if ( x == 0 )
				return 0;
			else if ( x > 0 )
				return exp(kappa_d[4] * x) - 1.0;
			else
				return exp(kappa_w[4] * x) - 1.0;
		}

		double cost_deriv_other(double d, double a) const {
			double x = scale_factor(a) * d;

			if ( x == 0 )
				return 0;
			else if ( x > 0 )
				return kappa_d[0] + pow(x / kappa_d[1], kappa_d[2]);
			else
				return -kappa_w[0] - pow(-x / kappa_w[1], kappa_w[2]);
		}

	public:
		AdjustmentCosts() = default;

		AdjustmentCosts(AdjustCostFnRatioMode mode_, bool exponential_costs_,
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
			}
			else {
				cost = [this](double d, double a) {return cost_fn_other(d, a);};
				cost1 = [this](double d, double a) {return cost_deriv_other(d, a);};
			}
		}
		
		std::function<double(double, double)> cost;
		std::function<double(double, double)> cost1;

		double cost1inv(double chi, double a) {return 0.0;}

		double scale_factor(double a) const {
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

		AdjustCostFnRatioMode mode;
		bool exponential_costs;
		double kappa_w_fc, kappa_d_fc;
		std::array<double, 5> kappa_w;
		std::array<double, 5> kappa_d;
};

#endif