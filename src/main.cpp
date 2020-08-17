
#include <iostream>
#include <procedures.h>
#include <parameters.h>
#include <options.h>
#include <model.h>
#include <steady_state.h>
#include <bellman.h>
#include <hank.h>
#include <Eigen/Dense>

int main () {
	std::string income_dir = "2point_3_5";

	Options options;
	Parameters params;

	Model model = Model(params, income_dir);
	
	SteadyState iss(model);

	HJB hjb(model, iss);

	// auto arr = new_array<double, 2>({3, 2});

	// for (int i=0; i<3; ++i)
	// 	for (int j=0; j<2; ++j)
	// 		arr[i][j] = i + j;

	std::cout << hjb.V[2][1][8] << '\n';
}
