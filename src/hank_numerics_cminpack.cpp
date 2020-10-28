#include <hank_numerics_cminpack.h>

#include <iostream>
#include <cminpack.h>
#include <cminpackP.h>

namespace {
	void check_cminpack_success(int info);
}

namespace HankNumerics {

void cminpack_hybrd1_wrapper(hank_cminpack_hybrd1_func_type fcn, void* args, int n, hank_float_type *x)
{
	hank_float_type fvec[n];
	double tol = 1.0e-9;

	int lwa = n * (3 * n + 13);
	hank_float_type wa[lwa];

	int info = cminpack_hybrd1_fnname(fcn, args, n, x, fvec, tol, wa, lwa);
	check_cminpack_success(info);
}

}

namespace {
	void check_cminpack_success(int info)
	{
		std::cout << '\n';
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
		std::cout << '\n';
	}
}
