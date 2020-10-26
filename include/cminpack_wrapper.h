#ifndef _CMINPACK_WRAPPER_H
#define _CMINPACK_WRAPPER_H

#include <hank_config.h>

#include <cminpack.h>
#include <cminpackP.h>

#include <utilities.h>

class Parameters;

class Model;

class Equilibrium;

class DistributionStatistics;

class IRF;

inline void check_cminpack_success(int info)
{
	std::cout << '\n';
	HankUtilities::horzline();
	HankUtilities::horzline();
	HankUtilities::horzline();
	if ( info == 0 ) {
		std::cout << "improper hybrd1 input parameters\n";
		throw 0;
	}
	else if ( info == 1 ) {
		std::cout << "hybrd1 has converged\n";
	}
	else if ( info == 2 ) {
		std::cout << "hybrd1 number of fcn calls has reached maximum\n";
		throw 0;
	}
	else if ( info == 3 ) {
		std::cout << "hybrd1 tol is too small, no further improvement possible\n";
		throw 0;
	}
	else if ( info == 4 ) {
		std::cout << "hybrd1 not making good progress\n";
		throw 0;
	}
	HankUtilities::horzline();
	HankUtilities::horzline();
	HankUtilities::horzline();
	std::cout << '\n';
}


template<typename T>
void cminpack_hybrd1_wrapper(cminpack_func_nn fcn, T* args, int n, real *x)
{
	real fvec[n];
	double tol = 1.0e-9;

	int lwa = n * (3 * n + 13);
	real wa[lwa];

	int info = cminpack_hybrd1_fnname(fcn, (void *) args, n, x, fvec, tol, wa, lwa);
	check_cminpack_success(info);
}

// template<typename F>
// class CminpackObjectiveFcn {
// 	public:
// 		CminpackObjectiveFcn(F& obj_fn_) : obj_fn(obj_fn_) {}

// 		F obj_fn;
// };

// template<typename T>
// void cminpack_hybrd1_wrapper_alt(T& objective, int n, real *x) {
// 	real fvec[n];
// 	double tol = 1.0e-9;

// 	int lwa = n * (3 * n + 13);
// 	real wa[lwa];

// 	int info = cminpack_hybrd1_fnname(hybrd1_obj_fn, (void *) objective, n, x, fvec, tol, wa, lwa);
// 	check_cminpack_success(info);
// }

// inline int hybrd1_obj_fn(void* fn_ptr, int /* n */, const real *x, real *fvec, int /* iflag */ ) {
// 	Hybrd1ObjectiveFn& obj_fn = *(Hybrd1ObjectiveFn *) fn_ptr;
// 	obj_fn(x, fvec);
// }

#endif