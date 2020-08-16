
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

	Options options;
	Parameters params;

	Model model = Model(params, income_dir);
	
	InitialSteadyState iss(model, options);
	iss = find_initial_steady_state(model);

	std::cout << model.get_rb_effective();
}
