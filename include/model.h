#ifndef _MODEL_H
#define _MODEL_H

#include <iostream>

#include <parameters.h>
#include <procedures.h>
#include <hank.h>

class Model {
	public:
		Model(Parameters params);

		vector bgrid;
		vector agrid;
		vector cgrid;
		vector occgrid;
		vector occdist;
		Parameters p;

		void make_grids();
};

#endif