
#include <iostream>
#include <procedures.h>
#include <parameters.h>
#include <model.h>
#include <initial_steady_state.h>
#include <hank.h>
#include <Eigen/Dense>

int main () {
	std::string income_dir = "2point_3_5";

	Parameters params = Parameters();
	Model model = Model(params, income_dir);
	
	solve_initial_steady_state(model);

}
