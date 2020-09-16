#include <stationary_dist.h>
#include <hank_eigen_dense.h>
#include <hank_eigen_sparse.h>
#include <model.h>
#include <steady_state.h>
#include <bellman.h>
#include <parameters.h>
#include <transition_matrix.h>
#include <iostream>
#include <assert.h>
#include <cmath>

#include <hank_macros.h>

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
	double_matrix lmat = speye(model.ny) + delta * model.ymarkovoff.transpose();
	int iabx = TO_INDEX_1D(0, p.nb_neg, p.na);

	double_matrix gmat = make_dist_guess(model, abdelta);
	double_matrix gmat_update(p.na * p.nb, model.ny);

	std::vector<sparse_matrix> B(model.ny);
	std::vector<sparse_solver> spsolvers(model.ny);
	for (int iy=0; iy<model.ny; ++iy) {
		sparse_matrix A = get_kfe_transition_matrix(p, model, ss.ra,
			hjb.optimal_decisions, iy);

		B[iy] = A.transpose();

		// Adjust A' matrix for non-linearly spaced grids
		B[iy] = inv_abdelta.asDiagonal() * B[iy] * abdelta.asDiagonal();
		B[iy] *= -delta;
		B[iy] += speye(p.na * p.nb) * (1.0 + delta * p.deathrate - delta * model.ymarkovdiag(iy,iy));
		B[iy].makeCompressed();

		spsolvers[iy].compute(B[iy]);
		if ( spsolvers[iy].info() != Eigen::Success )
				throw "Sparse solver failure";

	}

	double diff = 1.0e10;
	int ii = 0;
	while ( (diff > gtol) & (ii < maxiter) ) {
		for (int iy=0; iy<model.ny; ++iy) {
			double_vector lgmat = gmat * double_vector(lmat.row(iy));
			lgmat(iabx) += delta * p.deathrate * gmat.col(iy).dot(abdelta) / abdelta(iabx);

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
	gmat = gmat.array().max(0.0);

	double pmass = (gmat.array().colwise() * abdelta.array()).matrix().sum();
	assert( abs(1.0 - pmass) < 1.0e-6 );

	density = StdVector3d<double>(p.na, p.nb, model.ny);
	for (int ia=0; ia<p.na; ++ia)
		for (int ib=0; ib<p.nb; ++ib)
			for (int iy=0; iy<model.ny; ++iy)
				density(ia, ib, iy) = gmat(TO_INDEX_1D(ia, ib, p.na), iy);
}

namespace {
	double_matrix make_dist_guess(const Model& model, const double_vector& abdelta) {
		const Parameters& p = model.p;
		double gmass, p_y;
		double_matrix gmat = double_matrix::Zero(p.na * p.nb, model.ny);
		std::vector<int> vdims = {p.na, p.nb, model.ny};
		gmat.set_dims_3d(vdims.data(), 3);

		for (int iy=0; iy<model.ny; ++iy) {
			p_y = model.ydist(iy);
			if ( (p.deathrate == 0.0) & !p.borrowing ) {
				gmat.as3d(0, 1, iy) = p_y;
			}
			else if ( (p.deathrate == 0.0) & p.borrowing ) {
				gmat.as3d(0, p.nb_neg+1, iy) = p_y;
			}
			else if ( p.borrowing ) {
				gmat.as3d(0, p.nb_neg, iy) = p_y;
			}
			else {
				gmat.as3d(0, 0, iy) = p_y;
			}

			gmass = gmat.col(iy).dot(abdelta);
			gmat.col(iy) *= p_y / gmass;
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