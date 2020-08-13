#ifndef _MODEL_H
#define _MODEL_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>

#include <parameters.h>
#include <procedures.h>
#include <hank.h>

void read_vector(const std::string& file_loc, std::vector<double>& grid);

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

		Parameters p;

		void make_grids();
		void create_income_process(const std::string& income_dir);
};

#endif