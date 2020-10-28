#ifndef _HANK_H
#define _HANK_H

#include <hank_config.h>
#include <vector>
#include <iostream>
#include <memory>
#include <map>
#include <string>
#include <math.h>

#if STACK_LIQ_FIRST == 0
	#define TO_INDEX_1D_HANK_HEADER(a, b, y, na, nb) ((a) + (na) * (b) + (na) * (nb) * (y))
#elif STACK_LIQ_FIRST == 1
	#define TO_INDEX_1D_HANK_HEADER(a, b, y, na, nb) ((b) + (nb) * (a) + (nb) * (na) * (y))
#endif

namespace HANK {
	class OptimStatus
	{
		public:
			OptimStatus(const std::vector<std::string>& equation_names_,
				const std::vector<std::string>& variable_names_,
				const hank_float_type* deviations_,
				const std::vector<hank_float_type> variable_values_);

			void print() const;

			std::vector<std::string> equation_names, variable_names;

			std::vector<hank_float_type> deviations, variable_values;

			int n;
	};

	struct OptimNorm
	{
		OptimNorm(double norm_) : norm(norm_) {}

		double norm;
	};

	void print(const OptimStatus& optim_status);

	void print(const OptimNorm& optim_norm);

	const hank_float_type ValueNotSet = -91912395.1;

	void horzline();

	template<typename T>
	void print(const std::map<std::string, T>& variables, const std::string& title)
	{
		std::cout << '\n';
		horzline();
		std::cout << title << ":\n";
		for (auto variable : variables)
			std::cout << variable.first << " = " << variable.second << '\n';
		horzline();
		std::cout << '\n';
	}

	template<typename T>
	void print(const std::vector<T>& vec)
	{
		for (const auto& x : vec)
			std::cout << x << '\n';
	}

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

	template<typename T>
	T norm(const T* arr, int n)
	{
		T normval = 0;
		for (int i=0; i<n; ++i)
			normval += pow(arr[i], 2);

		normval = sqrt(normval);
		return normval;
	}
}

enum class AdjustCostFnRatioMode { none, linear, max };

template<typename T>
class StdVector3d {
	public:
		StdVector3d() {}

		StdVector3d(int n0, int n1, int n2) : shape{n0, n1, n2}, vector(n0 * n1 * n2) {}

		StdVector3d(const std::vector<int> dims_)
		{
			int n = 1;
			for (unsigned int i=0; i<dims_.size(); ++i) {
				shape[i] = dims_[i];

				if ( shape[i] > 0 )
					n *= shape[i];
			}

			vector.resize(n);
		}

		StdVector3d<T>& operator=(const StdVector3d<T>& other_vec)
		{
			for (int i=0; i<3; ++i)
				shape[i] = other_vec.shape[i];

			vector = other_vec.vector;

			return *this;
		}

		StdVector3d(const StdVector3d<T>& other_vec)
		{
			*this = other_vec;
		}

		operator std::vector<T>() const {return vector;}

		int shape[3];

		std::vector<T> vector;

		std::vector<T> as_vector() const {return vector;}

		T operator()(int i, int j, int k) const {return vector[TO_INDEX_1D_HANK_HEADER(i, j, k, shape[0], shape[1])];}

		T& operator()(int i, int j, int k) {return vector[TO_INDEX_1D_HANK_HEADER(i, j, k, shape[0], shape[1])];}

		int size() const {return vector.size();}

		T* data() {return vector.data();}

		const T* data() const {return vector.data();}
};

using vector3dr = StdVector3d<hank_float_type>;

#undef TO_INDEX_1D_HANK_HEADER

#endif