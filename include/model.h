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
		Parameters p;

		void make_grids();
};

#endif