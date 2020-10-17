
#include <hank_config.h>
#include <cminpack.h>
#include <cminpackP.h>

class Parameters;

class Model;

class Equilibrium;

class DistributionStatistics;

class IRF;

template<typename T>
int cminpack_hybrd1_wrapper(cminpack_func_nn fcn, T* args, int n, real *x) {
	real fvec[n];
	double tol = 1.0e-9;

	int lwa = n * (3 * n + 13);
	real wa[lwa];

	return cminpack_hybrd1_fnname(fcn, (void *) args, n, x, fvec, tol, wa, lwa);
}