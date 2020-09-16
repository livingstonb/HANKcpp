#include <distribution_statistics.h>
#include <parameters.h>
#include <model.h>
#include <bellman.h>
#include <stationary_dist.h>
#include <upwinding.h>
#include <iostream>
#include <algorithm>
#include <hank_numerics.h>
#include <utility>

#include <hank_macros.h>

namespace {
	void print_result(const std::string& expr, double val) {
		std::cout << expr << " = " << val << '\n';
	}

	VectorXi sort_by_values(VectorXd& vals, VectorXd& dist) {
		assert(vals.size() == dist.size());
		std::vector<std::pair<int, double>> zipped(vals.size());
		VectorXi indices(vals.size());
		for (unsigned int i=0; i<vals.size(); ++i)
			zipped[i] = std::make_pair(i, vals[i]);

		std::sort(zipped.begin(), zipped.end(),
			[](auto a, auto b) {return a.second < b.second;});

		VectorXd vcopy = vals;
		VectorXd dcopy = dist;

		for (unsigned int i=0; i<vals.size(); ++i)
			indices[i] = zipped[i].first;

		vals = vcopy(indices);
		dist = dcopy(indices);

		return indices;
	}
}

DistributionStatistics::DistributionStatistics(const Parameters& p_, const Model& model,
	const HJB& hjb, const StationaryDist& sdist) : p(p_) {

	const Upwinding::Policies& policies = hjb.optimal_decisions;

	ArrayXd abdelta(p.nab);

	MatrixXd nw_aby(p.nab, model.ny);
	MatrixXd agrid_aby(p.nab, model.ny);
	MatrixXd bgrid_aby(p.nab, model.ny);
	MatrixXd labor_aby(p.nab, model.ny);
	MatrixXd h_aby(p.nab, model.ny);

	int iab;
	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na);
			abdelta(iab) = model.adelta(ia) * model.bdelta(ib);
			for (int iy=0; iy<model.ny; ++iy) {
				nw_aby(iab, iy) = model.agrid(ia) + model.bgrid(ib);
				agrid_aby(iab, iy) = model.agrid(ia);
				bgrid_aby(iab, iy) = model.bgrid(ib);
				labor_aby(iab, iy) = policies.h(ia, ib, iy) * model.yprodgrid(iy);
				h_aby(iab, iy) = policies.h(ia, ib, iy);
			}
		}
	}

	// Distribution
	ArrayXd gdistvec = as_eigen<VectorXd>(sdist.density);
	MatrixXd gdistmat = Eigen::Map<MatrixXd>(gdistvec.data(), p.nab, model.ny);

	MatrixXd pdistmat = gdistmat.array().colwise() * abdelta;
	VectorXd pdistvec = eflatten(pdistmat);

	// Joint asset-income distributions
	MatrixXd p_ay = MatrixXd::Zero(p.na, model.ny);
	MatrixXd p_by = MatrixXd::Zero(p.nb, model.ny);

	iab = 0;
	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na);
			for (int iy=0; iy<model.ny; ++iy) {
				p_ay(ia, iy) += gdistmat(iab, iy) * abdelta(iab);
				p_by(ib, iy) += gdistmat(iab, iy) * abdelta(iab);
			}
		}
	}

	// Marginal distributions over illiq assets, liq assets
	VectorXd p_a = p_ay.rowwise().sum();
	VectorXd p_b = p_by.rowwise().sum();
	VectorXd pcum_a = cumsum(p_a);
	VectorXd pcum_b = cumsum(p_b);

	// Marginal distribution over net worth
	VectorXd nw_grid = nw_aby.col(0);
	VectorXd g_nw = gdistmat.rowwise().sum();
	VectorXi nw_order = sort_by_values(nw_grid, g_nw);

	VectorXd nwdelta(p.nab);
	nwdelta(seq(0, p.nab-2)) = 0.5 * (nw_grid(seq(1, p.nab-1)) - nw_grid(seq(0, p.nab-2)));
	nwdelta(seq(1, p.nab-1)) += 0.5 * (nw_grid(seq(1, p.nab-1)) - nw_grid(seq(0, p.nab-2)));

	VectorXd p_nw(p.nab);
	int inw;
	for (int iab=0; iab<p.nab; ++iab) {
		inw = nw_order(iab);
		p_nw(iab) = abdelta(inw) * gdistmat.row(inw).sum();
	}
	VectorXd pcum_nw = cumsum(p_nw);

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
		for (int ip=0; ip<model.nprod; ++ip) {
			Elabor_occ[io] += pdistmat.col(iy).dot(labor_aby.col(iy));
			++iy;
		}
	}

	// Percentiles
	double apct, bpct, nwpct;
	for (auto pval : pctiles) {
		// Illiquid asset
		if ( pcum_a(0) >= pval )
			apct = model.agrid[0];
		else
			apct = HankNumerics::lininterp1(
				p.na, pcum_a.data(), model.agrid.data(), pval);
		a_pctiles.push_back(apct);

		// Liquid asset
		if ( pcum_b(0) >= pval )
			bpct = model.bgrid[0];
		else
			bpct = HankNumerics::lininterp1(
				p.nb, pcum_b.data(), model.bgrid.data(), pval);
		b_pctiles.push_back(bpct);

		// Net worth
		if ( pcum_nw(0) >= pval )
			nwpct = nw_grid[0];
		else
			nwpct = HankNumerics::lininterp1(
				p.nab, pcum_nw.data(), nw_grid.data(), pval);
		nw_pctiles.push_back(nwpct);
	}
}


void DistributionStatistics::print() {
	std::cout << '\n' << '\n' << "--------------------------\n" << "Output:\n";
	print_result("E[h]", Ehours);
	print_result("E[nw]", Enetworth);
	print_result("E[a]", Ea);
	print_result("E[b]", Eb);
	print_result("Median(a)", a_pctiles[5]);
	print_result("Median(nw)", nw_pctiles[5]);
	print_result("Median(b)", b_pctiles[5]);
	std::cout << "--------------------------\n";
}
