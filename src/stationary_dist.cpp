#include <stationary_dist.h>
#include <hank_eigen_dense.h>
#include <hank_eigen_sparse.h>
#include <model.h>
#include <equilibrium.h>
#include <upwinding.h>
#include <parameters.h>
#include <transition_matrix.h>
#include <iostream>
#include <assert.h>
#include <math.h>

#include <hank_macros.h>

namespace {
	MatrixXr make_dist_guess(const Parameters& p, const Model& model);

	void check_progress(double vdiff, int freq, int ii, double vtol);

	void check_dist(const MatrixXr& distcheck, const Model& model);
}

void StationaryDist::compute(const Parameters& p, const Model& model, const Equilibrium& equm, const Upwinding::Policies& policies)
{
	if ( equm.is_transition_equilibrium() ) {
		dispfreq = 0;
		delta = equm.tdelta;
	}

	Eigen::Map<const VectorXr> abdeltavec(model.abdelta.data(), model.abdelta.size());
	VectorXr inv_abdelta = abdeltavec.cwiseInverse();
	int iabx = TO_INDEX_1D(0, p.nb_neg, p.na, p.nb);

	std::vector<SparseXd> B(model.ny);
	std::vector<sparse_solver> spsolvers(model.ny);

	if ( dispfreq > 0 )
		std::cout << "Beginning KFE iteration..." << '\n';

	for (int iy=0; iy<model.ny; ++iy) {
		SparseMatContainer Acont = get_kfe_transition_matrix(p, model, equm.ra, equm.illprice, equm.illpricedot,
			policies, iy);
		SparseXd& A = Acont.get();

		// Adjust A' matrix for non-linearly spaced grids
		B[iy] = inv_abdelta.cast<double>().asDiagonal() * A.transpose() * abdeltavec.cast<double>().asDiagonal();
		B[iy] *= -delta;
		B[iy] += speye(p.na * p.nb) * (1.0 + delta * p.deathrate - delta * model.matrices->ymarkovdiag(iy,iy));
		B[iy].makeCompressed();

		spsolvers[iy].compute(B[iy]);
		if ( spsolvers[iy].info() != Eigen::Success )
				throw "Sparse solver failure";

	}

	Eigen::MatrixXd lmat = deye(model.ny).cast<double>() + delta * model.matrices->ymarkovoff.cast<double>().transpose();
	MatrixXr gmat(p.nab, model.ny);
	MatrixXr gmat_update(p.nab, model.ny);
	if ( equm.is_transition_equilibrium() ) {
		for (int ia=0; ia<p.na; ++ia)
			for (int ib=0; ib<p.nb; ++ib)
				for (int iy=0; iy<model.ny; ++iy)
					gmat(TO_INDEX_1D(ia, ib, p.na, p.nb), iy) = density(ia, ib, iy);

		for (int iy=0; iy<model.ny; ++iy) {
			Eigen::VectorXd lgmat = gmat.cast<double>() * Eigen::VectorXd(lmat.row(iy));
			lgmat(iabx) += delta * p.deathrate * gmat.col(iy).dot(abdeltavec) / abdeltavec[iabx];

			gmat_update.col(iy) = spsolvers[iy].solve(lgmat).cast<hank_float_type>();
			if ( spsolvers[iy].info() != Eigen::Success )
				throw "Sparse solver failure";
		}

		gmat = gmat_update;
	}
	else {
		gmat = make_dist_guess(p, model);

		double diff = 1.0e10;
		int ii = 0;
		while ( (diff > gtol) & (ii < maxiter) ) {
			for (int iy=0; iy<model.ny; ++iy) {
				Eigen::VectorXd lgmat = gmat.cast<double>() * Eigen::VectorXd(lmat.row(iy));
				lgmat(iabx) += delta * p.deathrate * gmat.col(iy).dot(abdeltavec) / abdeltavec[iabx];

				gmat_update.col(iy) = spsolvers[iy].solve(lgmat).cast<hank_float_type>();
				if ( spsolvers[iy].info() != Eigen::Success )
					throw "Sparse solver failure";
			}

			diff = (gmat - gmat_update).cwiseAbs().maxCoeff();

			if ( dispfreq > 0 )
				check_progress(diff, dispfreq, ii, gtol);

			gmat = gmat_update;
			++ii;
		}

		if ( ii == maxiter ) {
			std::cout << "KFE did not converge" << '\n';
			throw 0;
		}
	}
	gmat = gmat.array().max(0.0);

	double pmass = (gmat.array().colwise() * abdeltavec.array()).matrix().sum();
	assert( abs(1.0 - pmass) < 1.0e-6 );

	density = vector3dr(p.na, p.nb, model.ny);
	for (int ia=0; ia<p.na; ++ia)
		for (int ib=0; ib<p.nb; ++ib)
			for (int iy=0; iy<model.ny; ++iy)
				density(ia, ib, iy) = gmat(TO_INDEX_1D(ia, ib, p.na, p.nb), iy);

	MatrixXr distcheck = gmat.array().colwise() * abdeltavec.array();
	check_dist(distcheck, model);
}

namespace {
	MatrixXr make_dist_guess(const Parameters& p, const Model& model)
	{
		MatrixXr gmat = MatrixXr::Zero(p.nab, model.ny);
		Eigen::Map<const VectorXr> abdeltavec(model.abdelta.data(), model.abdelta.size());

		double gmass;
		int bpos;
		for (int iy=0; iy<model.ny; ++iy) {
			if ( p.borrowing )
				bpos = p.nb_neg;
			else
				bpos = 0;

			if ( p.deathrate == 0.0 )
				++bpos;

			gmat(TO_INDEX_1D(0, bpos, p.na, p.nb), iy) = model.ydist[iy];
			gmass = gmat.col(iy).dot(abdeltavec);
			gmat.col(iy) *= model.ydist[iy] / gmass;
		}
		return gmat;
	}

	void check_progress(double gdiff, int freq, int ii, double gtol)
	{
		if ( ii == 0 )
			return;
		else if ( ii % freq == 0)  {
			std::cout << "Iteration " << ii << ", diff = " << gdiff << '\n';
		}

		if ( gdiff <= gtol )
			std::cout << "Converged after " << ii << " iterations." << '\n';
	}

	void check_dist(const MatrixXr& distcheck, const Model& model)
	{
		VectorXr py = distcheck.array().colwise().sum();
		for (int iy=0; iy<model.ny; ++iy) {
			assert( abs(py(iy) - model.ydist[iy]) < 1.0e-6 );
		}
	}
}