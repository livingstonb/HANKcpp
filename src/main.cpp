
#include <iostream>
#include <procedures.h>
#include <parameters.h>
#include <options.h>
#include <model.h>
#include <steady_state.h>
#include <hank.h>
#include <Eigen/Dense>

int main () {
	std::string income_dir = "2point_3_5";

	Parameters params;
	Model model = Model(params, income_dir);
	
	Options options;
	// InitialSteadyState iss(model, options);
	solve_initial_steady_state(model, options);

	std::cout << model.get_rb_effective();
}
