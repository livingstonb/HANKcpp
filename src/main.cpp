
#include <iostream>
#include <procedures.h>
#include <parameters.h>
#include <model.h>
#include <initial_steady_state.h>
#include <hank.h>
#include <Eigen/Dense>

int main () {
	Parameters params = Parameters();
	Model model = Model(params);
	
	solve_initial_steady_state(model);

	auto arr = new_array<3>({2,1,3});
	
	// array_shape<3> shape = {{2, 1, 3}};
	// array_type<3> arr(shape, boost::fortran_storage_order());

	arr[1][0][2] = 2.2;
	arr[1][0][1] = -1.5;

	for (int i=0; i<2; ++i)
		for (int j=0; j<1; ++j)
			for (int k=0; k<3; ++k)
				std::cout << arr[i][j][k] << '\n';

	array_type<3> arr2(arr[boost::indices[range()][range()][range(2,3)]]);
	auto map = to_eigen(arr2);
	std::cout << '\n' << map << '\n';
}
