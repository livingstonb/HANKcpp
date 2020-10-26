#ifndef _HANK_H
#define _HANK_H

#include <hank_config.h>
#include <vector>
#include <iostream>
#include <memory>
#include <map>
#include <string>

#if STACK_LIQ_FIRST == 0
	#define TO_INDEX_1D(a, b, y, na, nb) ((a) + (na) * (b) + (na) * (nb) * (y))
#elif STACK_LIQ_FIRST == 1
	#define TO_INDEX_1D(a, b, y, na, nb) ((b) + (nb) * (a) + (nb) * (na) * (y))
#endif

enum class AdjustCostFnRatioMode { none, linear, max };

namespace HANK {
	const hank_float_type ValueNotSet = -91912395.1;

	template<typename T>
	void print(const std::map<std::string, T>& variables, const std::string& title) {
		std::cout << '\n';
		horzline();
		std::cout << title << ":\n";
		for (auto variable : variables)
			std::cout << variable.first << " = " << variable.second << '\n';
		horzline();
		std::cout << '\n';
	}
}

using model_vector_type = std::vector<hank_float_type>;

class WealthTarget
{
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

template<typename T1, typename T2=void, typename T3=void, typename T4=void, typename T5=void>
class UniquePtrContainer {
	public:
		UniquePtrContainer() {}

		std::shared_ptr<T1> ptr1 = nullptr;

		std::shared_ptr<T2> ptr2 = nullptr;

		std::shared_ptr<T3> ptr3 = nullptr;

		std::shared_ptr<T4> ptr4 = nullptr;

		std::shared_ptr<T5> ptr5 = nullptr;
};

#undef TO_INDEX_1D

#endif