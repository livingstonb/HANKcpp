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
	double_matrix make_dist_guess(const Model& model);

	void check_progress(double vdiff, int freq, int ii, double vtol);
}

void StationaryDist::compute(const Model& model, const SteadyState& ss, const HJB& hjb) {
	const Parameters& p = model.p;

	double_vector inv_abdelta = model.abdelta.cwiseInverse();
	double_matrix lmat = speye(model.ny) + delta * model.ymarkovoff.transpose();
	int iabx = TO_INDEX_1D(0, p.nb_neg, p.na, p.nb);

	double_matrix gmat = make_dist_guess(model);
	double_matrix gmat_update(p.na * p.nb, model.ny);

	std::vector<sparse_solver> spsolvers(model.ny);
	for (int iy=0; iy<model.ny; ++iy) {
		sparse_matrix A = get_kfe_transition_matrix(p, model, ss.ra,
			hjb.optimal_decisions, iy);
		sparse_matrix B = A.transpose();

		// Adjust A' matrix for non-linearly spaced grids
		B = inv_abdelta.asDiagonal() * B * model.abdelta.asDiagonal();
		B *= -delta;
		B += (1.0 + delta * p.deathrate - delta * model.ymarkovdiag(iy,iy)) * speye(p.na * p.nb);
		B.makeCompressed();

		spsolvers[iy].compute(B);
		if ( spsolvers[iy].info() != Eigen::Success )
				throw "Sparse solver failure";

	}

	double diff = 1.0e10;
	int ii = 0;
	while ( (diff > gtol) & (ii < maxiter) ) {
		for (int iy=0; iy<model.ny; ++iy) {
			double_vector lgmat = gmat * double_vector(lmat.row(iy));
			lgmat(iabx) += delta * p.deathrate * gmat.col(iy).dot(model.abdelta) / model.abdelta(iabx);

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

	double pmass = (gmat.array().colwise() * model.abdelta.array()).matrix().sum();
	assert( abs(1.0 - pmass) < 1.0e-6 );

	density = StdVector3d<double>(p.na, p.nb, model.ny);
	for (int ia=0; ia<p.na; ++ia)
		for (int ib=0; ib<p.nb; ++ib)
			for (int iy=0; iy<model.ny; ++iy)
				density(ia, ib, iy) = gmat(TO_INDEX_1D(ia, ib, p.na, p.nb), iy);
}

namespace {
	MatrixXd make_dist_guess(const Model& model) {
		const Parameters& p = model.p;
		MatrixXd gmat = MatrixXd::Zero(p.na * p.nb, model.ny);

		double gmass, p_y;
		int bpos;
		for (int iy=0; iy<model.ny; ++iy) {
			p_y = model.ydist(iy);

			if ( p.borrowing )
				bpos = p.nb_neg;
			else
				bpos = 0;

			if ( p.deathrate == 0.0 )
				++bpos;

			gmat(TO_INDEX_1D(0, bpos, p.na, p.nb), iy) = p_y;
			gmass = gmat.col(iy).dot(model.abdelta);
			gmat.col(iy) *= p_y / gmass;
		}
		return gmat;
	}

	void check_progress(double gdiff, int freq, int ii, double gtol) {
		if ( ii == 0 )
			std::cout << "Beginning iteration" << '\n';
		else if ( ii % freq == 0)  {
			std::cout << "Iteration " << ii << ", diff = " << gdiff << '\n';
		}

		if ( gdiff <= gtol )
			std::cout << "Converged after " << ii << " iterations." << '\n';
	}
}