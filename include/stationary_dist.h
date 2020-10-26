#ifndef _STATIONARY_DIST_H
#define _STATIONARY_DIST_H

#include <vector>
#include <hank.h>

// Forward declarations
class Parameters;

class Model;

class HJB;

class Equilibrium;

class StationaryDist {
	public:
		StationaryDist() {}

		void compute(const Parameters& p, const Model& model, const Equilibrium& ss, const HJB& hjb);

		void transform();

		vector3dr density;

		int dispfreq = 50;

		double delta = 1.0e6;

		double gtol = 1.0e-8;

		int maxiter = 5000;
};

#endif