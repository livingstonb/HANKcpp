#ifndef _HANK_TYPES_H
#define _HANK_TYPES_H

#include <hank_config.h>
enum class LaborType { none, sep, ghh };

enum class AdjustCostFnRatioMode { none, linear, max };

enum class DepositCostMode { custom, symmetric, no_deposit_cost };

class Options {
	public:
		bool calibrateDiscountRate = false;
		bool equilibriumR = true;
		DepositCostMode depositCostMode = DepositCostMode::symmetric;
};

#endif