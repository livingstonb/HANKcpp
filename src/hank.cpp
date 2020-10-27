#include <hank.h>

namespace HANK {
	OptimStatus::OptimStatus(const std::vector<std::string>& equation_names_,
		const std::vector<std::string>& variable_names_,
		const hank_float_type* deviations_,
		const std::vector<hank_float_type> variable_values_)
	{
		equation_names = equation_names_;
		variable_names = variable_names_;
		variable_values = variable_values_;

		n = equation_names.size();
		for (int i=0; i<n; ++i)
			deviations.push_back(deviations_[i]);
	}

	void OptimStatus::print() const
	{
		std::cout << '\n';
		horzline();
		std::cout << "OPTIMIZATION RESULTS:\n";

		for (int i=0; i<n; ++i) {
			std::cout << variable_names[i] << " = " << variable_values[i];
			std::cout << ", ";
			std::cout << equation_names[i] << ": " << deviations[i];
			std::cout << '\n';
		}

		horzline();
		std::cout << '\n';
	}

	void print(const OptimStatus& optim_status)
	{
		optim_status.print();
	}

	void horzline()
	{
		std::cout << "--------------------------------\n";
	}
}