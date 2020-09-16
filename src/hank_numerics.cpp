#include <hank_numerics.h>
#include <math.h>

namespace HankNumerics {

double rtsec(std::function<double(double)> fn, double x1, double x2, double facc) {
	const int maxit = 20;
	double lrtsec, z1, z2;
	double dx, f, fl, xl, ltemp;

	z1 = x1;
	z2 = x2;
	fl = fn(z1);

	if (fabs(fl) < facc)
		return z1;

	f = fn(z2);

	if (fabs(fl) < fabs(f)) {
		lrtsec = z1;
		xl = z2;
		ltemp = fl;
		fl = f;
		f = ltemp;
	}
	else {
		xl = z1;
		lrtsec = z2;
	}

	for (int j=0; j<maxit; ++j) {
		dx = (xl - lrtsec) * f / (f - fl);
		xl = lrtsec;
		fl = f;
		lrtsec += dx;
		f = fn(lrtsec);

		if (fabs(f) < facc)
			return lrtsec;
	}

	throw "Exceeded maximum iterations";
}

double lininterp1(int n, const double *x, const double *y, double xi) {
	double xL, xH, yL, yH, maxel;
	int locL = -1;

	maxel = x[0];
	for (int i=1; i<n; ++i) {
		if ( (xi > x[i]) & (x[i] > maxel) ) {
			maxel = x[i];
			locL = i;
		}
	}

	if ( xi < x[0] )
		locL = 0;

	if ( locL >= n-1 )
		locL = n - 2;

	xL = x[locL];
	xH = x[locL+1];
	yL = y[locL];
	yH = y[locL+1];

	if ( abs(xL-xH < 1.0e-12) )
		return 0.5 * (yL + yH);
	else
		return yL  + ((xi - xL) / (xH - xL)) * (yH - yL);
}

}