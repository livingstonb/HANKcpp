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
	std::string markov_loc = "input/" + income_dir + "/ymarkov_combined.txt";
	read_matrix(grid_loc, logprodgrid);
	read_matrix(dist_loc, proddist);
	read_matrix(markov_loc, prodmarkov);

	for (auto x : logprodgrid)
		prodgrid.push_back(exp(x));
}

void read_matrix(const std::string& file_loc, vector& grid)
{
	std::string line, word;
	std::ifstream yfile;
	std::size_t current, previous;
	yfile.open(file_loc.data(), std::ios::in);



	while ( getline(yfile, line) ) {
		previous = 0;
    	current = find_multiple(line, 0);
    	if (current != std::string::npos) {
	    	while (current != std::string::npos) {
	    		if (current > 0) {
		    		word = line.substr(previous, current - previous);
		    		grid.push_back(std::stod(word));
			    }

			    previous = current + 1;
		    	current = find_multiple(line, previous);
		    }

		    word = line.substr(previous, current - previous);
		    grid.push_back(std::stod(word));
		}
		else {
			grid.push_back(std::stod(line));
		}
	}
	yfile.close();
}

std::size_t find_multiple(const std::string& line, int pos)
{
	std::size_t t1, t2;

	t1 = line.find("  ", pos);
	t2 = line.find(" -", pos);

	if ((t1 < t2) & (t1 != std::string::npos)) {
		return t1;
	}
	else {
		return t2;
	}
}