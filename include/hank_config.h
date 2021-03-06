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

// Precision (1 - double, 2 - long double)
#define HANK_PRECISION 2

// Precision for specific modules
// using model_float = double;

// --- DO NOT CHANGE ---
#if HANK_TURN_OFF_EIGEN_DEBUG == 1
	#define NDEBUG
	#define EIGEN_NO_STATIC_ASSERT
#endif

#define EIGEN_MATRIXBASE_PLUGIN "matrix_base_addons.h"

class Options {
	public:
		bool calibrateDiscountRate = false;
		bool equilibriumR = true;
		bool fast = false;
		bool print_diagnostics = false;
		bool skip_calibration = false;
};

extern const Options *global_hank_options;

#if HANK_PRECISION == 2
	using hank_float_type = long double;
	#define __cminpack_long_double__
	#define cminpack_hybrd1_fnname(args...) ldhybrd1(args)
#elif HANK_PRECISION == 1
	using hank_float_type = double;
	#define cminpack_hybrd1_fnname(args...) hybrd1(args)
#endif

#endif