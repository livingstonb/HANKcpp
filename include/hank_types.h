#ifndef _HANK_TYPES_H
#define _HANK_TYPES_H

#include <hank_config.h>
#include <vector>
#include <iostream>

#if STACK_LIQ_FIRST == 0
	#define TO_INDEX_1D(a, b, y, na, nb) ((a) + (na) * (b) + (na) * (nb) * (y))
#elif STACK_LIQ_FIRST == 1
	#define TO_INDEX_1D(a, b, y, na, nb) ((b) + (nb) * (a) + (nb) * (na) * (y))
#endif

enum class AdjustCostFnRatioMode { none, linear, max };

enum class DepositCostMode { custom, symmetric, no_deposit_cost };

class WealthTarget {
	public:
		enum class Type { mean, median, none };

		WealthTarget(WealthTarget::Type type_, double value_) : type(type_), value(value_) {}

		Type type;

		double value;

		bool is_mean() const {return (type == Type::mean);}

		bool is_median() const {return (type == Type::median);}
};

template<typename T>
class StdVector3d {
	public:
		StdVector3d() {}

		StdVector3d(int n0, int n1, int n2) : shape{n0, n1, n2}, vector(n0 * n1 * n2) {}

		StdVector3d(const std::vector<int> dims_) {
			int n = 1;
			for (unsigned int i=0; i<dims_.size(); ++i) {
				shape[i] = dims_[i];

				if ( shape[i] > 0 )
					n *= shape[i];
			}

			vector.resize(n);
		}

		operator std::vector<T>() const {return vector;}

		int shape[3];

		std::vector<T> vector;

		std::vector<T> as_vector() const {return vector;}

		T operator()(int i, int j, int k) const {
			return vector[TO_INDEX_1D(i, j, k, shape[0], shape[1])];
		}

		T& operator()(int i, int j, int k) {
			return vector[TO_INDEX_1D(i, j, k, shape[0], shape[1])];
		}

		int size() const {return vector.size();}

		T* data() {return vector.data();}

		const T* data() const {return vector.data();}
};

using vector3dr = StdVector3d<hank_float_type>;

#undef TO_INDEX_1D

#endif