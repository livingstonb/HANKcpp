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

std::vector<double> read_matrix(const std::string& file_loc);

std::size_t find_multiple(const std::string& line, int pos);

class Model
{
	public:
		Model(Parameters params, const std::string& income_dir);

		double_vector bgrid;
		double_vector agrid;
		double_vector cgrid;
		double_vector occgrid;
		double_vector occdist;

		double_vector logprodgrid;
		double_vector prodgrid;
		double_vector proddist;
		double_matrix prodmarkov;

		Parameters p;

		void make_grids();
		void create_income_process(const std::string& income_dir);
};

#endif