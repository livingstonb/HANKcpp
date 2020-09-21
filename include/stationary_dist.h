#ifndef _STATIONARY_DIST_H
#define _STATIONARY_DIST_H

#include <vector>
#include <hank_types.h>

// Forward declarations
class Model;

class HJB;

class SteadyState;

class StationaryDist {
	public:
		StationaryDist() {}

		void compute(const Model& model, const SteadyState& ss, const HJB& hjb);

		void transform();

		StdVector3d<double> density;

		int dispfreq = 50;

		double delta = 1.0e6;

		double gtol = 1.0e-8;

		int maxiter = 5000;
};

#endif