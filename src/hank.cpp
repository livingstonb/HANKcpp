#include <hank.h>

namespace {
	hank_float_type euc_norm(const hank_float_type* arr, int n);
}

std::map<std::string, hank_float_type> HankBase::variables_map() const
{
	return std::map<std::string, hank_float_type>();
}

std::string HankBase::title() const
{
	return "OBJECT VALUES";
}

namespace HANK {
	OptimStatus::OptimStatus(const std::vector<std::string>& equation_names_,
		const std::vector<std::string>& variable_names_,
		const hank_float_type* deviations_,
		const std::vector<hank_float_type>& variable_values_, int iter_)
	{
		equation_names = equation_names_;
		variable_names = variable_names_;
		variable_values = variable_values_;
		iter = iter_;

		n = equation_names.size();
		for (int i=0; i<n; ++i)
			deviations.push_back(deviations_[i]);

		norm = euc_norm(deviations_, n);
	}

	void OptimStatus::print() const
	{
		std::cout << "\titeration " << iter << "\n";

		for (int i=0; i<n; ++i) {
			std::cout << "\t  ";
			std::cout << variable_names[i] << " = " << variable_values[i];
			std::cout << ", ";
			std::cout << equation_names[i] << ": " << deviations[i];
			std::cout << "\n";
		}

		std::cout << "\t  norm = " << norm << '\n';
	}

	OptimNorm::OptimNorm(const hank_float_type* x, int n, int iter_)
		: iter(iter_)
	{
		norm = euc_norm(x, n);
	}

	void OptimNorm::print() const
	{
		std::cout << "\titeration " << iter << ", ";
		std::cout << "norm = " << norm << '\n';
	}

	void print(const HankBase* hank_obj)
	{	
		if ( hank_obj->override_print() )
			hank_obj->print();
		else
			print(hank_obj->variables_map(), hank_obj->title());
	}

	void horzline()
	{
		std::cout << "--------------------------------\n";
	}
}


namespace {
	hank_float_type euc_norm(const hank_float_type* arr, int n)
	{
		hank_float_type normval = 0;
		for (int i=0; i<n; ++i)
			normval += pow(arr[i], 2);
		
		normval = sqrt(normval);
		return normval;
	}
}
