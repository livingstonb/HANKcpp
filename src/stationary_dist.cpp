#include <stationary_dist.h>
#include <model.h>
#include <bellman.h>
#include <parameters.h>
#include <iostream>
#include <assert.h>
#include <cmath>

#define TO_INDEX_1D(a, b, na) ((a) + (na) * (b))

namespace {
	double_matrix make_dist_guess(const Model& model, const double_vector& abdelta);

	void check_progress(double vdiff, int freq, int ii, double vtol);
}

void StationaryDist::compute(const Model& model, const SteadyState& ss, const HJB& hjb) {
	const Parameters& p = model.p;
	double_vector abdelta(p.na * p.nb);
	for (int ia=0; ia<p.na; ++ia)
		for (int ib=0; ib<p.nb; ++ib)
			abdelta(TO_INDEX_1D(ia, ib, p.na)) = model.adelta(ia) * model.bdelta(ib);

	double_vector inv_abdelta = abdelta.cwiseInverse();
	double_matrix lmat = speye(p.ny) + delta * model.ymarkovoff.transpose();
	int iabx = TO_INDEX_1D(0, p.nb_neg, p.na);

	double_matrix gmat = make_dist_guess(model, abdelta);
	double_matrix gmat_update(p.na * p.nb, p.ny);

	std::vector<sparse_matrix> B(p.ny);
	std::vector<sparse_solver> spsolvers(p.ny);
	for (int iy=0; iy<p.ny; ++iy) {
		sparse_matrix A = hjb.get_A_matrix_KFE(ss, iy);

		B[iy] = A.transpose();

		// Adjust A' matrix for non-linearly spaced grids
		B[iy] = inv_abdelta.asDiagonal() * B[iy];
		B[iy] = -delta * B[iy] * abdelta.asDiagonal();
		B[iy] = B[iy] + speye(p.na * p.nb) * (1.0 + delta * p.deathrate - delta * model.ymarkovdiag(iy,iy));
		B[iy].makeCompressed();

		spsolvers[iy].compute(B[iy]);
		if ( spsolvers[iy].info() != Eigen::Success )
				throw "Sparse solver failure";

	}

	double diff = 1.0e10;
	int ii = 0;
	while ( (diff > gtol) & (ii < maxiter) ) {
		for (int iy=0; iy<p.ny; ++iy) {
			double_vector yvec_temp = lmat(iy,Eigen::all);
			map_type_vec yvec(yvec_temp.data(), yvec_temp.size());
			double_vector lgmat = gmat * yvec;
			lgmat(iabx) = lgmat(iabx) + delta * p.deathrate * gmat(Eigen::all,iy).dot(abdelta) / abdelta(iabx);

			gmat_update.col(iy) = spsolvers[iy].solve(lgmat);
			if ( spsolvers[iy].info() != Eigen::Success )
				throw "Sparse solver failure";
		}

		diff = (gmat - gmat_update).cwiseAbs().maxCoeff();
		check_progress(diff, dispfreq, ii, gtol);
		gmat = gmat_update;
		++ii;
	}

	if ( ii == maxiter )
		std::cout << "KFE did not converge" << '\n';

	double pmass = (gmat.array().colwise() * abdelta.array()).matrix().sum();
	assert(abs(1.0 - pmass) < 1.0e-6);

	density = gmat;
}

	

namespace {
	double_matrix make_dist_guess(const Model& model, const double_vector& abdelta) {
		const Parameters& p = model.p;
		double gmass, p_y;
		double_matrix gmat = double_matrix::Zero(p.na * p.nb, p.ny);
		std::vector<int> vdims = {p.na, p.nb, p.ny};
		gmat.set_dims_3d(vdims.data(), 3);

		for (int iy=0; iy<p.ny; ++iy) {
			p_y = model.ydist(iy);
			if ( (p.deathrate == 0.0) & !p.borrowing ) {
				// gmat.as3d(0, 1, iy) = p_y;
				gmat(TO_INDEX_1D(0, 1, p.na), iy) = p_y;
			}
			else if ( (p.deathrate == 0.0) & p.borrowing ) {
				// gmat.as3d(0, p.nb_neg+1, iy) = p_y;
				// gmat.as3d(1, p.nb_neg+1, iy) = p_y;
				gmat(TO_INDEX_1D(0, p.nb_neg+1, p.na), iy) = p_y;
				gmat(TO_INDEX_1D(1, p.nb_neg+1, p.na), iy) = p_y;
			}
			else if ( p.borrowing ) {
				// gmat.as3d(0, p.nb_neg, iy) = p_y;
				gmat(TO_INDEX_1D(0, p.nb_neg, p.na), iy) = p_y;
			}
			else {
				gmat(0, iy) = p_y;
			}

			gmass = (gmat.col(iy).cwiseProduct(abdelta)).sum();
			gmat.col(iy) = p_y * gmat.col(iy) / gmass;
		}

		return gmat;
	}

	void check_progress(double gdiff, int freq, int ii, double gtol) {
		if ( ii == 0 )
			std::cout << "Iteration " << ii << ", diff = " << gdiff << '\n';
		else if ( (ii == 1) | (ii % freq == 0) ) {
			std::cout << "Iteration " << ii << ", diff = " << gdiff << '\n';
		}

		if ( gdiff <= gtol )
			std::cout << "Converged after " << ii << " iterations." << '\n';
	}
}