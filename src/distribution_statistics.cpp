#include <distribution_statistics.h>
#include <parameters.h>
#include <model.h>
#include <bellman.h>
#include <stationary_dist.h>
#include <upwinding.h>
#include <iostream>
#include <algorithm>
#include <hank_numerics.h>

#define TO_INDEX_1D(a, b, na) ((a) + (na) * (b))

namespace {
	void print_result(const std::string& expr, double val) {
		std::cout << expr << " = " << val << '\n';
	}
}

DistributionStatistics::DistributionStatistics(const Parameters& p_, const Model& model,
	const HJB& hjb, const StationaryDist& sdist) : p(p_) {

	const Upwinding::Policies& policies = hjb.optimal_decisions;

	int nz = p.na * p.nb;
	ArrayXd abdelta(nz);

	MatrixXd nw_aby(p.na * p.nb, p.ny);
	MatrixXd agrid_aby(p.na * p.nb, p.ny);
	MatrixXd bgrid_aby(p.na * p.nb, p.ny);
	MatrixXd labor_aby(p.na * p.nb, p.ny);
	MatrixXd h_aby(p.na * p.nb, p.ny);

	int iab;
	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na);
			abdelta(iab) = model.adelta(ia) * model.bdelta(ib);
			for (int iy=0; iy<p.ny; ++iy) {
				nw_aby(iab, iy) = model.agrid(ia) + model.bgrid(ib);
				agrid_aby(iab, iy) = model.agrid(ia);
				bgrid_aby(iab, iy) = model.bgrid(ib);
				labor_aby(iab, iy) = policies.h(ia, ib, iy) * model.yprodgrid(iy);
				h_aby(iab, iy) = policies.h(ia, ib, iy);
			}
		}
	}

	ArrayXd gdistvec = as_eigen<VectorXd>(sdist.density);
	MatrixXd gdistmat = Eigen::Map<MatrixXd>(gdistvec.data(), p.na * p.nb, p.ny);

	MatrixXd pdistmat = gdistmat.array().colwise() * abdelta;
	VectorXd pdistvec = eflatten(pdistmat);

	MatrixXd p_ay = MatrixXd::Zero(p.na, p.ny);
	MatrixXd p_by = MatrixXd::Zero(p.nb, p.ny);

	iab = 0;
	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na);
			for (int iy=0; iy<p.ny; ++iy) {
				p_ay(ia, iy) += gdistmat(iab, iy) * model.bdelta(ib);
				p_by(ib, iy) += gdistmat(iab, iy) * model.adelta(ia);
			}
		}
	}

	VectorXd p_a = p_ay.rowwise().sum();
	VectorXd p_b = p_by.rowwise().sum();
	VectorXd pcum_a = cumsum(p_a);
	VectorXd pcum_b = cumsum(p_b);

	// Moments
	Ehours = pdistvec.dot(to_eigenv(policies.h));
	Enetworth = pdistvec.dot(eflatten(nw_aby));
	Ea = pdistvec.dot(eflatten(agrid_aby));
	Eb = pdistvec.dot(eflatten(bgrid_aby));

	Elabor_occ.resize(p.nocc);
	std::fill(Elabor_occ.begin(), Elabor_occ.end(), 0.0);

	int iy = 0;
	VectorXd labor_ab, pmass_ab;
	for (int io=0; io<p.nocc; ++io) {
		for (int ip=0; ip<p.nprod; ++ip) {
			Elabor_occ[io] += pdistmat.col(iy).dot(labor_aby.col(iy));
			++iy;
		}
	}

	// Percentiles
	for (auto pval : pctiles) {
		if ( pcum_a(0) >= pval )
			a_pctiles.push_back(model.agrid[0]);
		else {
			double apct = HankNumerics::lininterp1(
				p.na, pcum_a.data(), model.agrid.data(), pval);
			a_pctiles.push_back(apct);
		}
	}
}


void DistributionStatistics::print() {
	std::cout << '\n' << '\n' << "--------------------------\n" << "Output:\n";
	print_result("E[h]", Ehours);
	print_result("E[nw]", Enetworth);
	print_result("E[a]", Ea);
	print_result("E[b]", Eb);
	print_result("Median(a)", a_pctiles[5]);
	std::cout << "--------------------------\n";
}
