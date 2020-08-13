#include <model.h>

void ModelBase::make_grids(const Parameters& p) {
	bgrid_ = double_vector(nb_);
	agrid_ = double_vector(na_);
	powerSpacedGrid(nb_, p.bmin, p.bmax, p.bcurv, bgrid_);
	powerSpacedGrid(na_, p.amin, p.amax, p.acurv, agrid_);

	auto occgrids = occupationGrid(p);
	occgrid_ = vector2eigenv(occgrids.first);
	occdist_ = vector2eigenv(occgrids.second);
}

void ModelBase::create_income_process(const std::string& income_dir) {
	std::string grid_loc = "input/" + income_dir + "/ygrid_combined.txt";
	logprodgrid_ = vector2eigenv(read_matrix(grid_loc));

	std::string dist_loc = "input/" + income_dir + "/ydist_combined.txt";
	proddist_ = vector2eigenv(read_matrix(dist_loc));

	std::string markov_loc = "input/" + income_dir + "/ymarkov_combined.txt";
	int k = proddist_.size();
	prodmarkov_ = vector2eigenm(read_matrix(markov_loc), k, k);

	prodgrid_ = logprodgrid_.array().exp();
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