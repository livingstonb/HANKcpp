#include <hank_numerics.h>
#include <math.h>
#include <iostream>

#define HANK_INCLUDE_EIGEN_LU
#include <hank_eigen_dense.h>

namespace HankNumerics {

double rtsec(std::function<double(double)> fn, double x1, double x2, double facc)
{
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

double lininterp1(int n, const hank_float_type *x, const hank_float_type *y, double xi)
{
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

void jacobian_square(const broyden_fn_type& fn, int n, const hank_float_type *x,
	hank_float_type *f, hank_float_type *fjac, double step)
{
	for (int i=0; i<n; ++i)
		for (int j=0; j<n; ++j)
			fjac[i + n * j] = 0;

	hank_float_type xforjac[n];
	hank_float_type f1[n];
	for (int ix=0; ix<n; ++ix) {

		for (int i=0; i<n; ++i)
			xforjac[i] = x[i];

		xforjac[ix] += step;
		fn(n, xforjac, f1);

		for (int i=0; i<n; ++i)
			fjac[i + n * ix] = (f1[i] - f[i]) / step;
	}

	for (int i=0; i<n; ++i) {
		for (int j=0; j<n; ++j) {
			if ( fabs(fjac[i + n * j]) < 1.0e-7 )
				fjac[i + n * j] = 0;
		}
	}
}

void broyden_backstep(const broyden_fn_type& fn, int n, hank_float_type* x,
	hank_float_type* fvec, hank_float_type* fjac, int maxit, double ftol)
{
	bool recomputeJacobian = true;
	int nbs, maxbs = 5;
	double alpha = 1.0e-4;
	int lastitforjac = 0;
	double lstepmin = 1.0e-7;
	// int errorflag = 0;

	ArrayXr fvec0 = Eigen::Map<ArrayXr>(fvec, n);
	Eigen::Map<ArrayXr> fvecmap(fvec, n);
	ArrayXr x0 = Eigen::Map<ArrayXr>(x, n);
	Eigen::Map<ArrayXr> xmap(x, n);
	Eigen::Map<MatrixXr> fjacmap(fjac, n, n);

	double fnorm = 0.5 * fvec0.pow(2).sum();
	double fmaxerr = fvec0.abs().maxCoeff();

	if ( fmaxerr <= ftol )
		return;

	// MatrixXr fjacinv = invert_matrix(fjacmap, n, errorflag);
	MatrixXr fjacinv = fjacmap.inverse();

	ArrayXr ld, lgradf;
	MatrixXr ltemp1(n, 1);
	MatrixXr ltemp2(1, n);
	double fnorm0, fmaxerr0, lstep, lgradfld;
	for (int it=0; it<maxit; ++it) {
		// Newton direction
		ld = -fjacinv * fvecmap.matrix();

		// Store existing progress
		fvec0 = fvecmap;
		fnorm0 = fnorm;
		fmaxerr0 = fmaxerr;
		x0 = xmap;

		// Try full newton step
		lstep = 1.0;
	 	xmap = x0 + lstep * ld;
		fn(n, x, fvec);
		fnorm = 0.5 * fvecmap.pow(2).sum();
		fmaxerr = fvecmap.abs().maxCoeff();
		
		// First back step
		nbs = 1;
		// fjacmap = invert_matrix(fjacinv, n, errorflag);
		fjacmap = fjacinv.inverse();
		lgradf = fjacmap.transpose() * fvec0.matrix();
		lgradfld = lgradf.matrix().dot(lstep * ld.matrix());
		if ( fnorm > fnorm0 + alpha * lgradfld ) {
			++nbs;
			lstep = lgradfld / (2.0 * (fnorm - fnorm0 - lgradfld));
			xmap = x0 + lstep * ld;
			fn(n, x, fvec);
			fnorm = 0.5 * fvecmap.pow(2).sum();
			fmaxerr = fvecmap.abs().maxCoeff();
		}

		// Subsequent backsteps: use half steps
		for (int ibs=1; ibs<maxbs; ++ibs) {
			if ( fnorm > fnorm0 + alpha + lgradfld ) {
				++nbs;
				if ( nbs < maxbs )
					lstep *= 0.5;
				else if ( nbs == maxbs )
					lstep = fmin(lstepmin, 0.1 * lstep);

				xmap = x0 + lstep * ld;
				fn(n, x, fvec);
				fnorm = 0.5 * fvecmap.pow(2).sum();
				fmaxerr = fvecmap.abs().maxCoeff();
			}
		}

		if ( (nbs == maxbs) & (fnorm >= fnorm0) ) {
			xmap = x0;
			fvecmap = fvec0;

			if ( !recomputeJacobian ) {
				std::cerr << "Broyden: not a descent direction\n";
				return;
			}
			else {
				if ( it == lastitforjac ) {
					if ( it == 0 )
						std::cerr << "Broyden: not a descent direction\n";
					else if ( it > 1 )
						std::cerr << "Broyden: backstep did not work even after re-computing Jacobian\n";

					std::cerr << "Final:  L2 norm: " << sqrt(fnorm0) << ", max err: " << fmaxerr0 << "\n";
				}
				std::cout << "Broyden: recomputing Jacobian\n";
				jacobian_square(fn, n, x0.data(), fvec0.data(), fjac, 1.0e-8);
				// fjacinv = invert_matrix(fjac, n);
				fjacinv = fjacmap.inverse();
				lastitforjac = it;
			}
		}
		else {
			// Update estimate of inverse jacobian
			ArrayXr lu = fjacinv * (fvecmap - fvec0).matrix();
			ltemp1.col(0) = lstep * ld - lu;
			ltemp2.row(0) = (lstep * ld) * fjacinv.array();
			fjacinv += ltemp1 * ltemp2;
		}
	}
}


// MatrixXr invert_matrix(const MatrixXr& matrix, int n, int& errorflag)
// {
// 	MatrixXr inverse(n, n);
// 	MatrixXr augmatrix(n, 2 * n);
// 	bool flag = true;
// 	double m;

// 	// Augment input matrix with an identity matrix
// 	for (int i=0; i<n; ++i) {
// 		for (int j=0; j<2*n; ++j) {
// 			if ( j < n )
// 				augmatrix(i, j) = matrix(i, j);
// 			else if ( i + n == j )
// 				augmatrix(i, j) = 1;
// 			else
// 				augmatrix(i, j) = 0;
// 		}
// 	}

// 	// Reduce augmented matrix to upper traingular form
// 	for (int k=0; k<n-1; ++k) {
// 		if ( augmatrix(k,k) == 0 ) {
// 			flag = false;

// 			for (int i=k+1; i<n; ++i) {
// 				if ( augmatrix(i,k) != 0 ) {
// 					for (int j=0; j<2*n; ++j)
// 						augmatrix(k,j) += augmatrix(i,j);

// 					flag = true;
// 					break;
// 				}

// 				if ( !flag ) {
// 					inverse = MatrixXr::Zero(n,n);
// 					errorflag = -1;
// 					std::cerr << "Matrix is non-invertible\n";
// 					return inverse;
// 				}
// 			}
// 		}

// 		for (int j=k+1; j<n; ++j) {
// 			m = augmatrix(j,k) / augmatrix(k,k);

// 			for (int i=k; i<2*n; ++i) {
// 				augmatrix(j,i) = augmatrix(j,i) - m * augmatrix(k,i)
// 			}
// 		}
// 	}

// 	// Test for invertibility
// 	for (int i=0; i<n; ++i) {
// 		if ( augmatrix(i,i) == 0 ) {
// 			inverse = MatrixXr::Zero(n,n);
// 			errorflag = -1;
// 			std::cerr << "Matrix is non-invertible\n";
// 			return inverse;
// 		}
// 	}

// 	// Make diagonal elements as 1
// 	for (int i=0; i<n; ++i) {
// 		m = augmatrix(i,i);

// 		for (int j=i; j<2*n; ++j)
// 			augmatrix(i,j) /= m;
// 	}

// 	// Reduced right side half of augmented matrix to identity matrix
// 	for (int k=n-2; k>=0; --k) {
// 		for (int i=0; i<k; ++i) {
// 			m = augmatrix(i,k+1);
// 			for (int j=k; k<2*n; ++k)
// 				augmatrix(i,j) -= augmatrix(k+1,j) * m;
// 		}
// 	}

// 	// Store answer
// 	for (int i=0; i<n; ++i)
// 		for (int j=0; j<n; ++j)
// 			inverse(i,j) = augmatrix(i,j+n);

// 	errorflag = 0;
// 	return inverse;
// }

}