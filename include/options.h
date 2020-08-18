#ifndef _OPTIONS_H
#define _OPTIONS_H

#include <iostream>
#include <hank.h>

class Options {
	public:
		bool calibrateDiscountRate = false;
		bool equilibriumR = true;
		DepositCostMode depositCostMode = DepositCostMode::symmetric;
};

#endif