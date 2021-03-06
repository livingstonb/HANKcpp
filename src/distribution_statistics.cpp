#include <distribution_statistics.h>
#include <parameters.h>
#include <model.h>
#include <upwinding.h>
#include <iostream>
#include <algorithm>
#include <hank_numerics.h>

#include <hank_macros.h>
#include <hank_eigen_dense.h>

namespace {
	VectorXi sort_by_values(VectorXr& vals, VectorXr& dist);
}

DistributionStatistics::DistributionStatistics(const Parameters& p_, const Model& model,
	const Upwinding::Policies& policies, const vector3dr& density_) : density(density_)
{
	const Parameters& p = p_;

	MatrixXr nw_aby(p.nab, model.ny);
	MatrixXr agrid_aby(p.nab, model.ny);
	MatrixXr bgrid_aby(p.nab, model.ny);
	MatrixXr labor_aby(p.nab, model.ny);
	MatrixXr h_aby(p.nab, model.ny);

	int iab;
	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na, p.nb);
			for (int iy=0; iy<model.ny; ++iy) {
				nw_aby(iab, iy) = model.agrid[ia] + model.bgrid[ib];
				agrid_aby(iab, iy) = model.agrid[ia];
				bgrid_aby(iab, iy) = model.bgrid[ib];
				labor_aby(iab, iy) = policies.h(ia, ib, iy) * model.yprodgrid[iy];
				h_aby(iab, iy) = policies.h(ia, ib, iy);
			}
		}
	}

	// Distribution
	ArrayXr gdistvec = as_eigen<ArrayXr>(density);
	MatrixXr gdistmat = Eigen::Map<MatrixXr>(gdistvec.data(), p.nab, model.ny);

	VectorXr abdeltavec = as_eigen<VectorXr>(model.abdelta);
	MatrixXr pdistmat = gdistmat.array().colwise() * abdeltavec.array();
	VectorXr pdistvec = eflatten(pdistmat);

	// Joint asset-income distributions
	MatrixXr p_ay = MatrixXr::Zero(p.na, model.ny);
	MatrixXr p_by = MatrixXr::Zero(p.nb, model.ny);

	for (int ia=0; ia<p.na; ++ia) {
		for (int ib=0; ib<p.nb; ++ib) {
			iab = TO_INDEX_1D(ia, ib, p.na, p.nb);
			for (int iy=0; iy<model.ny; ++iy) {
				p_ay(ia, iy) += gdistmat(iab, iy) * model.abdelta[iab];
				p_by(ib, iy) += gdistmat(iab, iy) * model.abdelta[iab];
			}
		}
	}

	// Marginal distributions over illiq assets, liq assets
	VectorXr p_a = p_ay.rowwise().sum();
	VectorXr p_b = p_by.rowwise().sum();
	VectorXr pcum_a = cumsum(p_a);
	VectorXr pcum_b = cumsum(p_b);

	// Marginal distribution over net worth
	VectorXr nw_grid = nw_aby.col(0);
	VectorXr g_nw = gdistmat.rowwise().sum();
	VectorXi nw_order = sort_by_values(nw_grid, g_nw);

	VectorXr nwdelta(p.nab);
	nwdelta(seq(0, p.nab-2)) = 0.5 * (nw_grid(seq(1, p.nab-1)) - nw_grid(seq(0, p.nab-2)));
	nwdelta(seq(1, p.nab-1)) += 0.5 * (nw_grid(seq(1, p.nab-1)) - nw_grid(seq(0, p.nab-2)));

	VectorXr p_nw(p.nab);
	int inw;
	for (int iab=0; iab<p.nab; ++iab) {
		inw = nw_order(iab);
		p_nw(iab) = model.abdelta[inw] * gdistmat.row(inw).sum();
	}
	VectorXr pcum_nw = cumsum(p_nw);

	// Moments
	Ehours = pdistvec.dot(to_eigenv(policies.h));
	Enetworth = pdistvec.dot(eflatten(nw_aby));
	Ea = pdistvec.dot(eflatten(agrid_aby));
	Eb = pdistvec.dot(eflatten(bgrid_aby));
	Elabor = pdistvec.dot(eflatten(labor_aby));

	Elabor_occ.resize(p.nocc);
	std::fill(Elabor_occ.begin(), Elabor_occ.end(), 0.0);

	int iy = 0;
	double pocc;
	for (int io=0; io<p.nocc; ++io) {
		pocc = 0;
		for (int ip=0; ip<model.nprod; ++ip) {
			Elabor_occ[io] += pdistmat.col(iy).dot(labor_aby.col(iy));
			pocc += pdistmat.col(iy).sum();
			++iy;
		}
		Elabor_occ[io] /= pocc;
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

std::map<std::string, hank_float_type> DistributionStatistics::variables_map() const
{
	std::map<std::string, hank_float_type> variables;

	variables.insert({"E[h]", Ehours});
	variables.insert({"E[labor]", Elabor});
	variables.insert({"E[nw]", Enetworth});
	variables.insert({"E[a]", Ea});
	variables.insert({"E[b]", Eb});
	variables.insert({"Median(a)", a_pctiles[5]});
	variables.insert({"Median(b)", b_pctiles[5]});
	variables.insert({"Median(nw)", nw_pctiles[5]});

	return variables;
}

namespace {
	VectorXi sort_by_values(VectorXr& vals, VectorXr& dist)
	{
		assert(vals.size() == dist.size());
		std::vector<std::pair<int, double>> zipped(vals.size());
		VectorXi indices(vals.size());
		for (unsigned int i=0; i<vals.size(); ++i)
			zipped[i] = std::make_pair(i, vals[i]);

		std::sort(zipped.begin(), zipped.end(),
			[](auto a, auto b) {return a.second < b.second;});

		VectorXr vcopy = vals;
		VectorXr dcopy = dist;

		for (unsigned int i=0; i<vals.size(); ++i)
			indices[i] = zipped[i].first;

		vals = vcopy(indices);
		dist = dcopy(indices);

		return indices;
	}
}