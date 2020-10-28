#ifndef _HANK_NUMERICS_CMINPACK_H
#define _HANK_NUMERICS_CMINPACK_H

#include <hank_config.h>

typedef int(*hank_cminpack_hybrd1_func_type)(void *p, int n, const hank_float_type* x, hank_float_type* fvec, int iflag);

namespace HankNumerics {

void cminpack_hybrd1_wrapper(hank_cminpack_hybrd1_func_type fcn, void* args, int n, hank_float_type *x);
	
}

#endif