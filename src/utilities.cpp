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

}