
#include <hank_config.h>
#include <cminpack.h>
#include <cminpackP.h>

class Parameters;

class Model;

class EquilibriumElement;

class DistributionStatistics;

class IRF;

template<typename T>
int cmsinpack_hybrd1_wrapper(cminpack_func_nn fcn, T* args, int n, real *x) {
	real fvec[n];
	double tol = 1.0e-9;

	int lwa = n * (3 * n + 13);
	real wa[lwa];

	return cminpack_hybrd1_fnname(fcn, (void *) args, n, x, fvec, tol, wa, lwa);
}

template<typename T1, typename T2=void, typename T3=void, typename T4=void, typename T5=void>
class SolverArgs {
	public:
		SolverArgs(const T1* arg1_, int n)
			: arg1(arg1_) {x.reset(new hank_float_type[n]);}

		SolverArgs(const T1* arg1_, const T2* arg2_)
			: arg1(arg1_), arg2(arg2_) {x.reset(new hank_float_type[n]);}

		SolverArgs(const T1* arg1_, const T2* arg2_, const T3* arg3_)
			: arg1(arg1_), arg2(arg2_), arg3(arg3_) {x.reset(new hank_float_type[n]);}

		SolverArgs(const T1* arg1_, const T2* arg2_, const T3* arg3_, const T4* arg4_)
			: arg1(arg1_), arg2(arg2_), arg3(arg3_), arg4(arg4_) {x.reset(new hank_float_type[n]);}

		SolverArgs(const T1* arg1_, const T2* arg2_, const T3* arg3_, const T4* arg4_, const T5* arg5_)
			: arg1(arg1_), arg2(arg2_), arg3(arg3_), arg4(arg4_), arg5(arg5_) {x.reset(new hank_float_type[n]);}

		const T1* arg1 = nullptr;

		const T2* arg2 = nullptr;

		const T3* arg3 = nullptr;

		const T4* arg4 = nullptr;

		const T5* arg5 = nullptr;

		std::unique_ptr<hank_float_type[]> x = nullptr;
};