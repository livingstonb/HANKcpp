#include <distribution_statistics.h>
#include <parameters.h>
#include <model.h>
#include <bellman.h>
#include <stationary_dist.h>
#include <upwinding.h>
#include <iostream>
#include <algorithm>

#define TO_INDEX_1D(a, b, na) ((a) + (na) * (b))

namespace {
	void print_result(const std::string& expr, double val) {
		std::cout << expr << " = " << val << '\n';
	}
}

class DistributionStatistics::DistStruct {
	public:
		DistStruct(int nz, int ny)
			: networth(nz, ny), pmass(nz, ny),
				agrid_ab(nz, ny), bgrid_ab(nz, ny),
				wage(nz, ny), labor(nz, ny), pmass1d(nz) {}

		double_matrix networth, pmass, agrid_ab, bgrid_ab, wage, labor;

		double_vector pmass1d;
};

namespace {
	double take_expectation(const double_matrix& arr1, const double_matrix& arr2) {
		return (arr1.array() * arr2.array()).sum();
	}

	double take_expectation(const double_vector& arr1, const double_matrix& arr2) {
		return arr1.dot(eflatten(arr2));
	}
}

DistributionStatistics::DistributionStatistics(const Parameters& p_, const Model& model,
	const HJB& hjb, const StationaryDist& sdist) : p(p_) {

	int iab, nz = p.na * p.nb;
	DistStruct dist_struct(nz, p.ny);
	double_vector abdelta(nz);
	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na);
			abdelta(iab) = model.adelta(ia) * model.bdelta(ib);
			for (int iy=0; iy<p.ny; ++iy) {
				dist_struct.networth(iab, iy) = model.agrid(ia) + model.bgrid(ib);
				dist_struct.agrid_ab(iab, iy) = model.agrid(ia);
				dist_struct.bgrid_ab(iab, iy) = model.bgrid(ib);
				dist_struct.labor(iab, iy) = hjb.optimal_decisions.h(ia, ib, iy) * model.yprodgrid(iy);
				// dist_struct.wage(iab, iy) = model.yprodgrid(iy) * 
			}
		}
	}

	std::vector<double> density_copy = sdist.density.as_vector();
	dist_struct.pmass = map_type(density_copy.data(), p.na * p.nb, p.ny).array().colwise() * abdelta.array();
	dist_struct.pmass1d = map_type_vec(dist_struct.pmass.data(), p.na * p.nb * p.ny);
	compute_moments(hjb.optimal_decisions, dist_struct);
}

void DistributionStatistics::compute_moments(
	const Upwinding::Policies& policies, const DistStruct& dist_struct) {
	int iy;

	const double_matrix& pmass = dist_struct.pmass;
	auto expectation = [pmass](const auto& values) {
		return take_expectation(values, pmass);
	};

	Ehours = expectation(to_eigenv(policies.h));
	Enetworth = expectation(dist_struct.networth);
	Ea = expectation(dist_struct.agrid_ab);
	Eb = expectation(dist_struct.bgrid_ab);

	Elabor_occ.resize(p.nocc);
	std::fill(Elabor_occ.begin(), Elabor_occ.end(), 0.0);

	iy = 0;
	VectorXd labor_ab, pmass_ab;
	for (int io=0; io<p.nocc; ++io) {
		for (int ip=0; ip<p.nprod; ++ip) {
			labor_ab = dist_struct.labor.col(iy);
			pmass_ab = dist_struct.pmass.col(iy);
			Elabor_occ[io] += labor_ab.dot(pmass_ab);
			++iy;
		}
	}
}

void DistributionStatistics::print() {
	std::cout << '\n' << '\n' << "--------------------------\n" << "Output:\n";
	print_result("E[h]", Ehours);
	print_result("E[nw]", Enetworth);
	print_result("E[a]", Ea);
	print_result("E[b]", Eb);
	std::cout << "--------------------------\n";
}
