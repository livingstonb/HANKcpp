#include <distribution_statistics.h>
#include <hank_boost_eigen_routines.h>
#include <parameters.h>
#include <model.h>
#include <bellman.h>
#include <stationary_dist.h>
#include <upwinding.h>

#define TO_INDEX_1D(a, b, na) ((a) + (na) * (b))

class DistGrids {
	public:
		DistGrids(int nz, int ny) : networth(nz, ny), pmass(nz, ny), pmass1d(nz) {}
		double_matrix networth;
		double_matrix pmass;
		double_vector pmass1d;
};

namespace {
	double expectation(const double_matrix& arr1, const double_matrix& arr2) {
		return (arr1.array() * arr2.array()).sum();
	}
}

DistributionStatistics::DistributionStatistics(const Parameters& p, const Model& model,
	const HJB& hjb, const StationaryDist& sdist) {

	int iab;
	DistGrids grids(p.na * p.nb, p.ny);
	double_vector abdelta(p.na * p.nb);
	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na);
			abdelta(iab) = model.adelta(ia) * model.bdelta(ib);

			for (int iy=0; iy<p.ny; ++iy) {
				grids.networth(iab, iy) =  model.agrid(ia) + model.bgrid(ib);
			}
		}
	}

	std::vector<double> density_copy = sdist.density;
	grids.pmass = map_type(density_copy.data(), p.na * p.nb, p.ny).array().colwise() * abdelta.array();
	grids.pmass1d = map_type_vec(grids.pmass.data(), p.na * p.nb * p.ny);
	// grids.pmass = double_matrix(vdensity.colwise() * abdelta.array());
	// grids.pmass1d = double_vector(eflatten(grids.pmass));

	compute_moments(model, hjb.optimal_decisions, grids);
}

void DistributionStatistics::compute_moments(const Model& model,
	const Upwinding::Policies& policies, const DistGrids& grids) {
	Ehours = boost_dot(policies.h, grids.pmass1d);
	Enetworth = expectation(grids.networth, grids.pmass);
}