#include <transition_matrix.h>
#include <parameters.h>
#include <model.h>
#include <bellman.h>
#include <upwinding.h>
#include <hank_eigen_dense.h>

#define TO_INDEX_1D(a, b, na) ((a) + (na) * (b))

sparse_matrix construct_transition_matrix(const Parameters& p, const Model& model, double ra,
	const Upwinding::Policies& policies, int iy, bool kfe) {

	double d, s, acost, areturn, val, val1, val2;
	int iab;
	Bellman::Drifts drifts;

	double_vector adriftvec = (ra + p.perfectAnnuityMarkets * p.deathrate) * model.agrid.array();

	int na = p.na;
	triplet_list Aentries;
	Aentries.reserve(5 * p.na * p.nb);

	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na);
			d = policies.d(ia,ib,iy);
			s = policies.s(ia,ib,iy);
			acost = model.adjcosts.cost(d, model.agrid(ia));
			areturn = adriftvec(ia);

			// Compute drifts
			drifts = Bellman::Drifts(s, d, areturn, acost, kfe);

			// Matrix entries, ia-1
			if ( (ia > 0) & (drifts.aB != 0.0) ) {
				val = -drifts.aB / model.dagrid(ia-1);
				Aentries.push_back(triplet_type(iab, TO_INDEX_1D(ia-1, ib, na), val));
			}

			// Matrix entries, ib-1
			if ( (ib > 0) & (drifts.bB != 0.0) ) {
				val = -drifts.bB / model.dbgrid(ib-1);
				Aentries.push_back(triplet_type(iab, TO_INDEX_1D(ia, ib-1, na), val));
			}

			// Matrix entries, diagonal
			if ( ia == 0 )
				val1 = -drifts.aF / model.dagrid(ia);
			else if ( ia == p.na - 1 )
				val1 = drifts.aB / model.dagrid(ia-1);
			else
				val1 = drifts.aB / model.dagrid(ia-1) - drifts.aF / model.dagrid(ia);

			if ( ib == 0 )
				val2 = -drifts.bF / model.dbgrid(ib);
			else if ( ib == p.nb - 1 )
				val2 = drifts.bB / model.dbgrid(ib-1);
			else
				val2 = drifts.bB / model.dbgrid(ib-1) - drifts.bF / model.dbgrid(ib);

			if ( val1 + val2 != 0 )
				Aentries.push_back(triplet_type(iab, iab, val1 + val2));

			// Matrix entries, ia+1
			if ( (ia < p.na - 1 ) & (drifts.aF != 0.0) ) {
				val = drifts.aF / model.dagrid(ia);
				Aentries.push_back(triplet_type(iab, TO_INDEX_1D(ia+1, ib, na), val));
			}
			
			// Matrix entries, ib+1
			if ( (ib < p.nb - 1) & (drifts.bF != 0.0) ) {
				val = drifts.bF /  model.dbgrid(ib);
				Aentries.push_back(triplet_type(iab, TO_INDEX_1D(ia, ib+1, na), val));
			}
		}
	}

	sparse_matrix A = sparse_matrix(p.na * p.nb, p.na * p.nb);
	A.setFromTriplets(Aentries.begin(), Aentries.end());
	return A;
}

sparse_matrix get_kfe_transition_matrix(const Parameters& p, const Model& model, double ra,
	const Upwinding::Policies& policies, int iy) {
	bool kfe = true;
	return construct_transition_matrix(p, model, ra, policies, iy, kfe);
}