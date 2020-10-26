#include <transition_matrix.h>
#include <parameters.h>
#include <model.h>
#include <bellman.h>
#include <upwinding.h>
#include <hank_eigen_dense.h>
#include <hank_eigen_sparse.h>
#include <hank_macros.h>
#include <adjustment_costs.h>

namespace {
	struct Drifts {
		Drifts() {}

		Drifts(double s, double d, double areturn, double acost, bool kfe, double illprice);

		double aB, aF, bB, bF;
	};
}

SparseMatContainer construct_transition_matrix(const Parameters& p, const Model& model, double ra,
	double illprice, double illpricedot, const Upwinding::Policies& policies, int iy, bool kfe)
{

	double d, s, acost, areturn, val, val1, val2;
	int iab;
	Drifts drifts;

	auto agridvec = as_eigen_map<const ArrayXr>(model.agrid);
	VectorXr adriftvec = (ra - illpricedot / illprice + p.perfectAnnuityMarkets * p.deathrate) * agridvec;

	std::vector<EigenTriplet> Aentries;
	Aentries.reserve(5 * p.na * p.nb);

	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na, p.nb);
			d = policies.d(ia,ib,iy);
			s = policies.s(ia,ib,iy);
			acost = model.adjcosts->cost(d, model.agrid[ia]);
			areturn = adriftvec(ia);

			// Compute drifts
			drifts = Drifts(s, d, areturn, acost, kfe, illprice);

			// Matrix entries, ia-1
			if ( (ia > 0) & (drifts.aB != 0.0) ) {
				val = -drifts.aB / model.dagrid[ia-1];
				Aentries.push_back(EigenTriplet(iab, TO_INDEX_1D(ia-1, ib, na, p.nb), val));
			}

			// Matrix entries, ib-1
			if ( (ib > 0) & (drifts.bB != 0.0) ) {
				val = -drifts.bB / model.dbgrid[ib-1];
				Aentries.push_back(EigenTriplet(iab, TO_INDEX_1D(ia, ib-1, na, p.nb), val));
			}

			// Matrix entries, diagonal
			if ( ia == 0 )
				val1 = -drifts.aF / model.dagrid[ia];
			else if ( ia == p.na - 1 )
				val1 = drifts.aB / model.dagrid[ia-1];
			else
				val1 = drifts.aB / model.dagrid[ia-1] - drifts.aF / model.dagrid[ia];

			if ( ib == 0 )
				val2 = -drifts.bF / model.dbgrid[ib];
			else if ( ib == p.nb - 1 )
				val2 = drifts.bB / model.dbgrid[ib-1];
			else
				val2 = drifts.bB / model.dbgrid[ib-1] - drifts.bF / model.dbgrid[ib];

			if ( val1 + val2 != 0 )
				Aentries.push_back(EigenTriplet(iab, iab, val1 + val2));

			// Matrix entries, ia+1
			if ( (ia < p.na - 1 ) & (drifts.aF != 0.0) ) {
				val = drifts.aF / model.dagrid[ia];
				Aentries.push_back(EigenTriplet(iab, TO_INDEX_1D(ia+1, ib, na, p.nb), val));
			}
			
			// Matrix entries, ib+1
			if ( (ib < p.nb - 1) & (drifts.bF != 0.0) ) {
				val = drifts.bF /  model.dbgrid[ib];
				Aentries.push_back(EigenTriplet(iab, TO_INDEX_1D(ia, ib+1, na, p.nb), val));
			}
		}
	}

	SparseXd A = SparseXd(p.na * p.nb, p.na * p.nb);
	A.setFromTriplets(Aentries.begin(), Aentries.end());
	return SparseMatContainer(std::move(A));
}

SparseMatContainer get_kfe_transition_matrix(const Parameters& p, const Model& model, double ra,
	double illprice, double illpricedot, const Upwinding::Policies& policies, int iy) {
	bool kfe = true;
	return construct_transition_matrix(p, model, ra, illprice, illpricedot, policies, iy, kfe);
}

namespace {
	Drifts::Drifts(double s, double d, double areturn, double acost, bool kfe, double illprice)
	{
		if ( kfe ) {
			aB = fmin(d / illprice + areturn, 0.0);
			aF = fmax(d / illprice + areturn, 0.0);
			bB = fmin(s - d - acost, 0.0);
			bF = fmax(s - d - acost, 0.0);
		}
		else {
			aB = fmin(d / illprice, 0.0) + fmin(areturn, 0.0);
			aF = fmax(d / illprice, 0.0) + fmax(areturn, 0.0);
			bB = fmin(-d - acost, 0) + fmin(s, 0.0);
			bF = fmax(-d - acost, 0) + fmax(s, 0.0);
		}
	}
}