#include <hank_numerics.h>
#include <math.h>
#include <iostream>
#include <utilities.h>

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

	if (fabs(f) < facc)
		return z2;

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

	std::cout << "Exceeded maximum iterations";
	throw 0;
}

double lininterp1(int n, const hank_float_type *x, const hank_float_type *y, double xi) {
	hank_float_type xL, xH, yL, yH, maxel;
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

void jacobian_square(const broyden_fn_type& fn, int n, const hank_float_type *x, hank_float_type *f, hank_float_type *fjac, double step) {
	HankUtilities::fillarr(fjac, 0.0, n, n);

	hank_float_type* xforjac = new hank_float_type[n];
	hank_float_type* f1 = new hank_float_type[n];
	for (int ix=0; ix<n; ++ix) {
		HankUtilities::copyarr(x, xforjac, n);
		xforjac[ix] += step;
		fn(n, xforjac, f1);

		for (int ii=0; ii<n; ++ii)
			fjac[ii + n * ix] = (f1[ii] - f[ii]) / step;
	}

	for (int i=0; i<n; ++i) {
		for (int j=0; j<n; ++j) {
			if ( fabs(fjac[i + n * j]) < 1.0e-7 )
				fjac[i + n * j] = 0;
		}
	}

	delete[] xforjac;
	delete[] f1;
}

}