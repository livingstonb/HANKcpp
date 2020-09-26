#ifndef _HANK_NUMERICS_H
#define _HANK_NUMERICS_H

#include <hank_config.h>
#include <functional>

namespace HankNumerics {

double rtsec(std::function<double(double)> fn, double x1, double x2, double facc);

template<typename T, typename V, typename R>
double lininterp1(int n, const T *x, const V *y, R xi) {
	double xL, xH, yL, yH, maxel;
	int locL = -1;

	maxel = -1e12;
	for (int i=1; i<n; ++i) {
		if ( (xi > x[i]) & (x[i] > maxel) ) {
			maxel = x[i];
			locL = i;
		}
	}

	if ( xi <= x[0] )
		locL = 0;

	if ( locL >= n-1 )
		locL = n - 2;

	xL = x[locL];
	xH = x[locL+1];
	yL = y[locL];
	yH = y[locL+1];

	if ( abs(xL-xH) < 1.0e-12 )
		return 0.5 * (yL + yH);
	else
		return yL  + ((xi - xL) / (xH - xL)) * (yH - yL);
}

}

#endif