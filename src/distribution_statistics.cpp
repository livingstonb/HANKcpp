#include <distribution_statistics.h>
#include <hank_boost_eigen_routines.h>
#include <parameters.h>
#include <model.h>
#include <bellman.h>
#include <stationary_dist.h>
#include <upwinding.h>

#define TO_INDEX_1D(a, b, na) ((a) + (na) * (b))

namespace {
	double expectation(const double_matrix& arr1, const double_matrix& arr2) {
		return (arr1.array() * arr2.array()).sum();
	}
}

DistributionStatistics::DistributionStatistics(const Parameters& p, const Model& model,
	const HJB& hjb, const StationaryDist& sdist) {

	int iab;
	grids.networth = double_matrix(p.na * p.nb, p.ny);
	double_vector abdelta(p.na * p.nb);
	for (int ia=0; ia<p.na; ++ia) {

		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na);
			abdelta(TO_INDEX_1D(ia, ib, p.na)) = model.adelta(ia) * model.bdelta(ib);

			for (int iy=0; iy<p.ny; ++iy) {
				grids.networth(iab, iy) = model.bgrid(ib) + model.agrid(ia);
			}
		}
	}

	pmass = sdist.density.array().colwise() * abdelta.array();
	pmass1d = eflatten(pmass);
	compute_moments(model, hjb.optimal_decisions, grids);
}

void DistributionStatistics::compute_moments(const Model& model,
	const Upwinding::Policies& policies, const DistGrids& grids) {
	Ehours = boost_dot(policies.h, pmass1d);
	Enetworth = expectation(grids.networth, pmass);
}