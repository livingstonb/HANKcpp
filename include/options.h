#ifndef _OPTIONS_H
#define _OPTIONS_H

#include <iostream>

enum class DepositCostMode { custom, symmetric, no_deposit_cost };

class Options {
	public:
		bool calibrateDiscountRate = false;
		bool equilibriumR = true;
		DepositCostMode depositCostMode = DepositCostMode::symmetric;
};

#endif