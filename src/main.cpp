
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

	boost_array_type<double, 2> arr(boost::extents[2][3]);

	for (int i=0; i<2; ++i) {
		for (int j=0; j<3; ++j) {
			arr[i][j] = i + j;
		}
	}

	auto earr = boost2eigen(arr);
	std::cout << earr << '\n';

	boost_array_type<double, 2> temp_arr = reshape_array(arr, {6, 1});

	// auto earr2 = boost2eigen(temp_arr);
	std::cout << boost2eigen(temp_arr) << '\n';


	hjb.iterate(iss);

}
