#ifndef _HANK_NUMERICS_H
#define _HANK_NUMERICS_H

#include <hank_config.h>
#include <functional>

namespace HankNumerics {

double rtsec(std::function<double(double)> fn, double x1, double x2, double facc);

double lininterp1(int n, const hank_float_type *x, const hank_float_type *y, double xi);

template<typename T>
void linspace(double x, double y, int n, T& cont) {
	for (int i=0; i<n; ++i)
		cont[i] = x + i * (y - x) / (n - 1);
}

using broyden_fn_type = std::function<void(int, const hank_float_type*, hank_float_type*)>;

void jacobian_square(const broyden_fn_type& fn, int n, const hank_float_type *x,
	const hank_float_type *z, hank_float_type *fjac, double step);

void broyden_backstep(const broyden_fn_type& fn, int n, hank_float_type* x,
	hank_float_type* fvec, hank_float_type* fjac, int maxit, double ftol);

// MatrixXr invert_matrix(const MatrixXr& matrix, int n, int& errorflag);

}

#endif