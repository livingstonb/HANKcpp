#ifndef _HANK_NUMERICS_H
#define _HANK_NUMERICS_H

#include <hank_config.h>
#include <functional>

namespace HankNumerics {

double rtsec(std::function<double(double)> fn, double x1, double x2, double facc);

double lininterp1(int n, const double *x, const double *y, double xi);

}

#endif