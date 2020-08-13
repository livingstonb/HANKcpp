#include <model.h>

Model::Model(Parameters params, const std::string& income_dir) {
	p = params;
	make_grids();
	create_income_process(income_dir);
}

void Model::make_grids() {
	bgrid = double_vector(p.nb);
	agrid = double_vector(p.na);
	powerSpacedGrid(p.nb, p.bmin, p.bmax, p.bcurv, bgrid);
	powerSpacedGrid(p.na, p.amin, p.amax, p.acurv, agrid);

	auto occgrids = occupationGrid(p);
	occgrid = vector2eigenv(occgrids.first);
	occdist = vector2eigenv(occgrids.second);
}

void Model::create_income_process(const std::string& income_dir) {
	std::string grid_loc = "input/" + income_dir + "/ygrid_combined.txt";
	logprodgrid = vector2eigenv(read_matrix(grid_loc));

	std::string dist_loc = "input/" + income_dir + "/ydist_combined.txt";
	proddist = vector2eigenv(read_matrix(dist_loc));

	std::string markov_loc = "input/" + income_dir + "/ymarkov_combined.txt";
	int k = proddist.size();
	prodmarkov = vector2eigenm(read_matrix(markov_loc), k, k);

	prodgrid = logprodgrid.array().exp();
}

std::vector<double> read_matrix(const std::string& file_loc)
{
	std::string line, word;
	std::ifstream yfile;
	std::size_t current, previous;
	yfile.open(file_loc.data(), std::ios::in);

	std::vector<double> out;

	while ( getline(yfile, line) ) {
		previous = 0;
    	current = find_multiple(line, 0);
    	if (current != std::string::npos) {
	    	while (current != std::string::npos) {
	    		if (current > 0) {
		    		word = line.substr(previous, current - previous);
		    		out.push_back(std::stod(word));
			    }

			    previous = current + 1;
		    	current = find_multiple(line, previous);
		    }

		    word = line.substr(previous, current - previous);
		    out.push_back(std::stod(word));
		}
		else {
			out.push_back(std::stod(line));
		}
	}
	yfile.close();

	return out;
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