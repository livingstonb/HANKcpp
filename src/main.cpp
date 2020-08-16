
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
	
	SteadyState iss(model);

	std::cout << iss.capital << '\n';
}
