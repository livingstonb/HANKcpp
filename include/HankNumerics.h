#ifndef _RTSEC_H
#define _RTSEC_H

#include<functional>
#include<math.h>

namespace HankNumerics {

double rtsec(std::function<double(double)> fn, double x1, double x2, double facc);

}

#endif