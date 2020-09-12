#include <stationary_dist.h>
#include <hank_eigen_dense.h>
#include <hank_eigen_sparse.h>
#include <model.h>
#include <bellman.h>
#include <parameters.h>
#include <iostream>

#define TO_INDEX_1D(a, b, na) ((a) + (na) * (b))

namespace {
	double_matrix make_dist_guess(const Model& model, const double_vector& abdelta);
}

void StationaryDist::compute(const Model& model, const HJB& hjb, const std::vector<sparse_matrix>& A) {
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
	for (int iy=0; iy<p.ny; ++iy) {
		sparse_matrix B = A[iy].transpose();

		// Adjust A' matrix for non-linearly spaced grids
		B = inv_abdelta.asDiagonal() * B;
		B = B * abdelta.asDiagonal();

		B = -delta * B;
		B = B + speye(p.ny) * (1.0 + delta * p.deathrate - delta * model.ymarkovdiag(iy,iy));
		B.makeCompressed();

		double_vector lgmat = gmat * lmat(iy,Eigen::all);
		lgmat(iabx) = lgmat(iabx) + delta * p.deathrate * gmat(Eigen::all,iy).dot(abdelta);

		sparse_solver solver;
		solver.compute(B);
		if ( solver.info() != Eigen::Success )
			throw "Sparse solver failure";

		gmat_update.col(iy) = solver.solve(lgmat);
		if ( solver.info() != Eigen::Success )
			throw "Sparse solver failure";
	}
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
			if ( (p.deathrate == 0.0) & p.borrowing )
				gmat.as3d(0, 1, iy) = p_y;
			else if ( (p.deathrate == 0.0) & !p.borrowing ) {
				gmat.as3d(0, p.nb_neg+1, iy) = p_y;
				gmat.as3d(1, p.nb_neg+1, iy) = p_y;
			}
			else if ( p.borrowing )
				gmat(0, iy) = p_y;
			else
				gmat.as3d(0, p.nb_neg, iy) = p_y;

			gmass = (gmat.col(iy).cwiseProduct(abdelta)).sum();
			gmat.col(iy) = p_y * gmat.col(iy) / gmass;
		}

		return gmat;
	}
}