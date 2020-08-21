
#include <iostream>
#include <parameters.h>
#include <hank_types.h>
#include <model.h>
#include <steady_state.h>
#include <bellman.h>
#include <adjustment_costs.h>

int main () {
	std::string income_dir = "2point_3_5";

	Options options;

	Parameters params;
	params.setup(options);

	Model model = Model(params, income_dir);

	SteadyState iss(model);

	HJB hjb(model, iss); 
	hjb.iterate(iss);


}
