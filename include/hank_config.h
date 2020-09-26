#ifndef _HANK_CONFIG_H
#define _HANK_CONFIG_H

//  ---CONFIG OPTIONS ---
// Turn off Eigen debugging, set to 1 for performance
#define HANK_TURN_OFF_EIGEN_DEBUG 0

// Choice of sparse linear solver
// 0 - QR
// 1 - LU
// 2 - UmfPackLU
#define HANK_EIGEN_SPARSE_SOLVER 2

// Put liquid asset in first dimension instead of illiquid asset, for stacked variables
#define STACK_LIQ_FIRST 1

// Precision (0 - float, 1 - double, 2 - long double)
#define HANK_PRECISION 2

// --- DO NOT CHANGE ---
#if HANK_TURN_OFF_EIGEN_DEBUG == 1
	#define NDEBUG
#endif

#define EIGEN_MATRIXBASE_PLUGIN "matrix_base_addons.h"

class Options {
	public:
		bool calibrateDiscountRate = false;
		bool equilibriumR = true;
		bool fast = false;
		bool print_diagnostics = false;
};

extern const Options *global_hank_options;

#if HANK_PRECISION == 2
	using fp_type = long double;
	#define __cminpack_long_double__
// 	#define HANK_CMINPACK_SOLVER ldhybrd1
// #elif HANK_PRECISION == 1
// 	#define HANK_CMINPACK_SOLVER hybrd1
#elif HANK_PRECISION == 1
	using fp_type = double;
#endif

#endif