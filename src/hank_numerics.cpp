#include <hank_numerics.h>
#include<math.h>

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
		lrtsec = lrtsec + dx;
		f = fn(lrtsec);

		if (fabs(f) < facc)
			return lrtsec;
	}

	throw "Exceeded maximum iterations";
}

}