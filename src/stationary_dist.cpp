#include <stationary_dist.h>
#include <hank_eigen_dense.h>
#include <model.h>
#include <bellman.h>
#include <parameters.h>
#include <iostream>

#define TO_INDEX_1D(a, b, na) ((a) + (na) * (b))

namespace {
	double_matrix make_dist_guess(const Model& model, const double_vector& abdelta);
}

void StationaryDist::compute(const Model& model, const HJB& hjb) {
	const Parameters& p = model.p;
	double_vector abdelta(p.na * p.nb);
	for (int ia=0; ia<p.na; ++ia)
		for (int ib=0; ib<p.nb; ++ib)
			abdelta(TO_INDEX_1D(ia, ib, p.na)) = model.adelta(ia) * model.bdelta(ib);

	double_matrix gmat = make_dist_guess(model, abdelta);
}

	

namespace {
	double_matrix make_dist_guess(const Model& model, const double_vector& abdelta) {
		const Parameters& p = model.p;
		double gmass, p_y;
		double_matrix gmat = double_matrix::Zero(p.na * p.nb, p.ny);
		gmat.set_dims_3d({p.na, p.nb, p.ny});

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