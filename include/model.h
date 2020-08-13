#ifndef _MODEL_H
#define _MODEL_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include <parameters.h>
#include <procedures.h>
#include <hank.h>

void read_matrix(const std::string& file_loc, std::vector<double>& grid);

std::size_t find_multiple(const std::string& line, int pos);

class Model
{
	public:
		Model(Parameters params, const std::string& income_dir);

		vector bgrid;
		vector agrid;
		vector cgrid;
		vector occgrid;
		vector occdist;

		vector logprodgrid;
		vector prodgrid;
		vector proddist;
		vector prodmarkov;

		Parameters p;

		void make_grids();
		void create_income_process(const std::string& income_dir);
};

#endif