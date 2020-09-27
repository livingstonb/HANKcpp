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
	MatrixXr make_dist_guess(const Model& model);

	void check_progress(double vdiff, int freq, int ii, double vtol);

	void check_dist(const MatrixXr& distcheck, const Model& model) {
		VectorXr py = distcheck.array().colwise().sum();
		for (int iy=0; iy<model.ny; ++iy) {
			assert( abs(py(iy) - model.ydist(iy)) < 1.0e-6 );
		}
	}
}

void StationaryDist::compute(const Model& model, const SteadyState& ss, const HJB& hjb) {
	const Parameters& p = model.p;

	Eigen::VectorXd inv_abdelta = model.abdelta.cast<double>().cwiseInverse();
	Eigen::MatrixXd lmat = deye(model.ny).cast<double>() + delta * model.ymarkovoff.cast<double>().transpose();
	int iabx = TO_INDEX_1D(0, p.nb_neg, p.na, p.nb);

	MatrixXr gmat = make_dist_guess(model);
	MatrixXr gmat_update(p.nab, model.ny);

	std::vector<sparse_matrix> B(model.ny);
	std::vector<sparse_solver> spsolvers(model.ny);
	for (int iy=0; iy<model.ny; ++iy) {
		SparseMatContainer Acont = get_kfe_transition_matrix(p, model, ss.ra,
			hjb.optimal_decisions, iy);
		SparseXd& A = Acont.matrix;

		// Adjust A' matrix for non-linearly spaced grids
		B[iy] = inv_abdelta.asDiagonal() * A.transpose() * model.abdelta.cast<double>().asDiagonal();
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
			Eigen::VectorXd lgmat = gmat.cast<double>() * Eigen::VectorXd(lmat.row(iy));
			lgmat(iabx) += delta * p.deathrate * gmat.col(iy).dot(model.abdelta) / model.abdelta(iabx);

			gmat_update.col(iy) = spsolvers[iy].solve(lgmat).cast<hank_float_type>();
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

	density = vector3dr(p.na, p.nb, model.ny);
	for (int ia=0; ia<p.na; ++ia)
		for (int ib=0; ib<p.nb; ++ib)
			for (int iy=0; iy<model.ny; ++iy)
				density(ia, ib, iy) = gmat(TO_INDEX_1D(ia, ib, p.na, p.nb), iy);

	MatrixXr distcheck = gmat.array().colwise() * model.abdelta.array();
	check_dist(distcheck, model);
}

namespace {
	MatrixXr make_dist_guess(const Model& model) {
		const Parameters& p = model.p;
		MatrixXr gmat = MatrixXr::Zero(p.nab, model.ny);

		double gmass;
		int bpos;
		for (int iy=0; iy<model.ny; ++iy) {
			if ( p.borrowing )
				bpos = p.nb_neg;
			else
				bpos = 0;

			if ( p.deathrate == 0.0 )
				++bpos;

			gmat(TO_INDEX_1D(0, bpos, p.na, p.nb), iy) = model.ydist(iy);
			gmass = gmat.col(iy).dot(model.abdelta);
			gmat.col(iy) *= model.ydist(iy) / gmass;
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