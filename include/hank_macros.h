#ifndef _HANK_MACROS_H
#define _HANK_MACROS_H

#include <hank_config.h>

#if STACK_LIQ_FIRST == 0
	#define TO_INDEX_1D(a, b, na, nb) ((a) + (na) * (b))
#elif STACK_LIQ_FIRST == 1
	#define TO_INDEX_1D(a, b, na, nb) ((b) + (nb) * (a))
#endif

#endif