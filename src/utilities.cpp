#include <utilities.h>

namespace HankUtilities {

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

void check_cminpack_success(int info) {
	std::cout << '\n';
	HankUtilities::horzline();
	HankUtilities::horzline();
	HankUtilities::horzline();
	if ( info == 0 ) {
		std::cout << "improper hybrd1 input parameters\n";
		throw 0;
	}
	else if ( info == 1 ) {
		std::cout << "hybrd1 has converged\n";
	}
	else if ( info == 2 ) {
		std::cout << "hybrd1 number of fcn calls has reached maximum\n";
		throw 0;
	}
	else if ( info == 3 ) {
		std::cout << "hybrd1 tol is too small, no further improvement possible\n";
		throw 0;
	}
	else if ( info == 4 ) {
		std::cout << "hybrd1 not making good progress\n";
		throw 0;
	}
	HankUtilities::horzline();
	HankUtilities::horzline();
	HankUtilities::horzline();
	std::cout << '\n';
}

}