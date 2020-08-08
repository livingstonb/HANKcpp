#ifndef _MODEL_H
#define _MODEL_H

#include <iostream>

#include <parameters.h>
#include <procedures.h>

class Model {
	public:
		Model(Parameters params);

		std::vector<double> bgrid;
		std::vector<double> agrid;
		std::vector<double> occgrid;
		std::vector<double> occdist;
		Parameters p;

		void make_grids();
};

#endif