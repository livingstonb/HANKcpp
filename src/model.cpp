#include <model.h>

Model::Model(Parameters params, const std::string& income_dir) {
	p = params;
	make_grids();
	create_income_process(income_dir);
}

void Model::make_grids() {
	bgrid = powerSpacedGrid<vector>(p.nb, p.bmin, p.bmax, p.bcurv);
	agrid = powerSpacedGrid<vector>(p.na, p.amin, p.amax, p.acurv);

	auto occgrids = occupationGrid<vector>(p);
	occgrid = occgrids.first;
	occdist = occgrids.second;
}

void Model::create_income_process(const std::string& income_dir) {
	std::string grid_loc = "input/" + income_dir + "/ygrid_combined.txt";
	std::string dist_loc = "input/" + income_dir + "/ydist_combined.txt";
	read_vector(grid_loc, logprodgrid);
	read_vector(dist_loc, proddist);

	for (auto x : logprodgrid)
		prodgrid.push_back(exp(x));
}

void read_vector(const std::string& file_loc, vector& grid)
{
	std::string line, whole_file = "";
	std::ifstream yfile;
	yfile.open(file_loc.data(), std::ios::in);

	while ( getline (yfile, line) ) {
		grid.push_back( std::stod(line) );
		whole_file += line + '\n';
	}
	yfile.close();
}
