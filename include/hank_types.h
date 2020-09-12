#ifndef _HANK_TYPES_H
#define _HANK_TYPES_H

#include <hank_config.h>
#include <vector>

enum class LaborType { none, sep, ghh };

enum class AdjustCostFnRatioMode { none, linear, max };

enum class DepositCostMode { custom, symmetric, no_deposit_cost };

class Options {
	public:
		bool calibrateDiscountRate = false;
		bool equilibriumR = true;
		DepositCostMode depositCostMode = DepositCostMode::symmetric;
		bool fast = false;
};

template<typename T>
class StdVector3d {
	public:
		StdVector3d() {}

		StdVector3d(int n0, int n1, int n2) : dims({n0, n1, n2}), vector(n0 * n1 * n2) {}

		int dims[3];

		std::vector<T> vector;

		T operator()(int i, int j, int k) const {
			return vector[i + dims[0] * j + dims[0] * dims[1] * k];
		}

		T& operator()(int i, int j, int k) {
			return vector[i + dims[0] * j + dims[0] * dims[1] * k];
		}

		int size() const {return vector.size();}

		T* data() {return vector.data();}

		const T* data() const {return vector.data();}
};

#endif