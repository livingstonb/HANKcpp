#ifndef _STATIONARY_DIST_H
#define _STATIONARY_DIST_H

#include <hank_eigen_dense.h>
#include <hank_eigen_sparse.h>

// Forward declarations
class Model;

class HJB;

class SteadyState;

class StationaryDist {
	public:
		StationaryDist() {}

		void compute(const Model& model, const SteadyState& ss, const HJB& hjb);

		double_matrix density;
		int dispfreq = 1;
		double delta = 1.0e6;
		double gtol = 1.0e-12;
		int maxiter = 2000;
};

#endif