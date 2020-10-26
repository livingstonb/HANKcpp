#include <hank.h>
#include <model.h>
#include <utilities.h>
#include <hank_numerics.h>
#include <functions.h>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <assert.h>
#include <hank_macros.h>
#include <hank_eigen_dense.h>
#include <parameters.h>
#include <adjustment_costs.h>

namespace {
	void make_asset_grids(Model* model, const Parameters& p);

	void make_occupation_grids(Model* model, const Parameters& p);

	void create_income_process(Model* model, const Parameters& p);

	void create_combined_variables(Model* model, const Parameters& p);

	void check_nbl(double bgrid0, const Parameters& p);

	void fix_rounding(MatrixXr& mat);

	std::vector<hank_float_type> compute_grid_deltas(
		const std::vector<hank_float_type>& grid, const std::vector<hank_float_type>& dgrid);

	void powerSpacedGrid(double low, double high, double curv, std::vector<hank_float_type>& grid);

	void adjustPowerSpacedGrid(std::vector<hank_float_type>& grid);

	void print_value(const std::string& pname, double value, bool insert_endline);

	void print_value(const std::string& pname, double value);

	void check_adjcosts(const Parameters& p, const std::shared_ptr<AdjustmentCosts>& adjcosts);
}

Model::Model(const Parameters& p_) : p(p_) {
	matrices.reset(new ModelMatrices);

	make_asset_grids(this, p);
	make_occupation_grids(this, p);
	create_income_process(this, p);
	create_combined_variables(this, p);

	nocc = p.nocc;
	naby = p.nb * p.na * nocc * nprod;
	ntot = p.nb * p.na * nocc * nprod;
	dims = std::vector<int>({p.na, p.nb, ny});

	assertions();

	check_nbl(bgrid[0], p);

	double adjcost1max = 1.0e30;

	AdjustmentCosts adjcosts_initial(p.adjCostRatioMode, p.exponential_adjcosts,
		p.kappa_w_fc, p.kappa_d_fc, p.kappa_w, p.kappa_d, adjcost1max, p.dmax);

	adjcost1max = adjcosts_initial.cost1(p.dmax, 1.0);
	adjcosts.reset(new AdjustmentCosts(p.adjCostRatioMode, p.exponential_adjcosts,
		p.kappa_w_fc, p.kappa_d_fc, p.kappa_w, p.kappa_d, adjcost1max, p.dmax));

	check_adjcosts(p, adjcosts);

	if ( global_hank_options->print_diagnostics )
		print_values();
}

std::vector<hank_float_type> Model::get_rb_effective(hank_float_type rb, hank_float_type rborr) const {
	std::vector<hank_float_type> rb_effective;
	rb_effective.reserve(p.nb);
	for (auto el : bgrid) {
		if ( el >= 0.0 )
			rb_effective.push_back(rb + p.perfectAnnuityMarkets * p.deathrate);
		else
			rb_effective.push_back(rborr + p.perfectAnnuityMarkets * p.deathrate);
	}

	return rb_effective;
}

hank_float_type Model::util(hank_float_type c, hank_float_type riskaver) const {
	return HankFunctions::utility(c, p.prefshock, riskaver);
}

hank_float_type Model::util1(hank_float_type c, hank_float_type riskaver) const {
	return HankFunctions::utility1(c, p.prefshock, riskaver);
}

hank_float_type Model::util1inv(hank_float_type u, hank_float_type riskaver) const {
	return HankFunctions::utility1inv(u, p.prefshock, riskaver);
}

hank_float_type Model::labdisutil(hank_float_type h, hank_float_type chi) const {
	return HankFunctions::labor_disutility(h, p.frisch, chi);
}

hank_float_type Model::labdisutil1(hank_float_type h, hank_float_type chi) const {
	return HankFunctions::labor_disutility1(h, p.frisch, chi);
}

hank_float_type Model::labdisutil1inv(hank_float_type du, hank_float_type chi) const {
	return HankFunctions::labor_disutility1inv(du, p.frisch, chi);
}

hank_float_type Model::capadjcost(hank_float_type x) const {
	return HankFunctions::capadjcost(x, p.capadjcost, p.depreciation);
}

hank_float_type Model::capadjcost1(hank_float_type x) const {
	return HankFunctions::capadjcost1(x, p.capadjcost, p.depreciation);
}


hank_float_type Model::capadjcost1inv(hank_float_type x) const {
	return HankFunctions::capadjcost1inv(x, p.capadjcost, p.depreciation);
}

double Model::util1BC(double h, double riskaver, double chi, double bdrift, double netwage, double wagescale) const {
	double c = bdrift + h * netwage;
	assert( c > 0 );
	return labdisutil1(h, chi) - util1(c, riskaver) * netwage * p.labwedge / wagescale;
}

void Model::print_values() const {
	HankUtilities::horzline();
	std::cout << "COMPUTED VALUES, MODEL OBJECT:\n";
	print_value("nocc", nocc);
	print_value("nprod", nprod);
	print_value("ny", ny);

	for (int io=0; io<nprod; ++io)
		print_value("prodgrid[io]", prodgrid[io]);

	for (int io=0; io<nprod; ++io)
		print_value("proddist[io]", proddist[io]);

	print_value("E[prod]", EigenFunctions::dot(prodgrid, proddist));
	
	HankUtilities::horzline();
}

void Model::assertions() const {
	if ( abs(1.0 - EigenFunctions::sum(proddist)) > 1.0e-7 ) {
		std::cerr << "Proddist does not sum to one\n";
		throw 0;
	}

	if ( abs(1.0 - EigenFunctions::sum(occdist)) > 1.0e-7 ) {
		std::cerr << "Occdist does not sum to one\n";
		throw 0;
	}

	if ( abs(1.0 - EigenFunctions::sum(ydist)) > 1.0e-7 ) {
		std::cerr << "ydist does not sum to one\n";
		throw 0;
	}

	ArrayXr rowsums = matrices->ymarkov.rowwise().sum();
	if ( rowsums.abs().maxCoeff() > 1.0e-7 ) {
		std::cerr << "Markov ytrans matrix rows do not sum to zero\n";
		throw 0;
	}
}

namespace {
	void make_asset_grids(Model* model, const Parameters& p) {
		// Liquid asset
		model->bgrid.resize(p.nb);
		if ( !p.borrowing ) {
			powerSpacedGrid(p.bmin, p.bmax, p.bcurv, model->bgrid);
		}
		else {
			std::vector<hank_float_type> bgridpos(p.nb_pos);
			powerSpacedGrid(0.0, p.bmax, p.bcurv, bgridpos);

			double nbl = -p.lumptransfer / (p.rborr + p.perfectAnnuityMarkets * p.deathrate);
			double abl = fmax(nbl + p.cmin, p.blim);

			int nn = static_cast<int>(floor(p.nb_neg / 2.0) + 1);
			std::vector<hank_float_type> bgridneg(nn);
			powerSpacedGrid(abl, (abl+bgridpos[0])/2.0, p.bcurv_neg, bgridneg);

			std::copy(bgridneg.begin(), bgridneg.end(), model->bgrid.begin());
			std::copy(bgridpos.begin(), bgridpos.end(), model->bgrid.begin() + p.nb_neg);

			for (int i=nn; i<p.nb_neg; ++i)
				model->bgrid[i] = model->bgrid[p.nb_neg] - (model->bgrid[p.nb_neg-i] - model->bgrid[0]);
		}
		
		for (int ib=0; ib<p.nb-1; ++ib)
			model->dbgrid.push_back(model->bgrid[ib+1] - model->bgrid[ib]);

		model->bdelta = compute_grid_deltas(model->bgrid, model->dbgrid);

		// Illiquid asset
		model->agrid.resize(p.na);
		if ( p.oneAssetNoCapital ) {
			std::fill(model->agrid.begin(), model->agrid.end(), 0);
			std::fill(model->dagrid.begin(), model->dagrid.end(), 1);
			std::fill(model->adelta.begin(), model->adelta.end(), 1);
		}
		else {
			powerSpacedGrid(p.amin, p.amax, p.acurv, model->agrid);
			adjustPowerSpacedGrid(model->agrid);

			for (int ia=0; ia<p.na-1; ++ia)
				model->dagrid.push_back(model->agrid[ia+1] - model->agrid[ia]);

			model->adelta = compute_grid_deltas(model->agrid, model->dagrid);
		}

		model->abdelta.resize(p.nab);
		for (int ia=0; ia<p.na; ++ia)
			for (int ib=0; ib<p.nb; ++ib)
				model->abdelta[TO_INDEX_1D(ia, ib, p.na, p.nb)] = model->adelta[ia] * model->bdelta[ib];
	}

	void make_occupation_grids(Model* model, const Parameters& p) {
		// Occupation types
		if ( p.nocc == 1 ) {
			model->occYsharegrid = std::vector<hank_float_type>({1});
			model->occNsharegrid = std::vector<hank_float_type>({1});
			model->occdist = std::vector<hank_float_type>({1});
		}
		else if ( p.nocc == 4 ) {
			model->occYsharegrid = std::vector<hank_float_type>({0.325, 0.275, 0.225, 0.175});
			model->occNsharegrid = std::vector<hank_float_type>({0.1, 0.15, 0.25, 0.5});
			model->occdist = std::vector<hank_float_type>({0.25, 0.25, 0.25, 0.25});
		}
		else {
			std::cerr << "Invalid number of occupation points\n";
			throw 0;
		}
		model->nocc = p.nocc;
	}

	void create_income_process(Model* model, const Parameters& p) {
		std::string grid_loc = "../input/" + p.income_dir + "/ygrid_combined.txt";
		model->logprodgrid = HankUtilities::read_matrix(grid_loc);

		std::string dist_loc = "../input/" + p.income_dir + "/ydist_combined.txt";
		model->proddist = HankUtilities::read_matrix(dist_loc);

		std::string markov_loc = "../input/" + p.income_dir + "/ymarkov_combined.txt";

		if ( p.adjustProdGridFrisch )
			for (auto& el : model->logprodgrid)
				el /= (1.0 + p.adjFrischGridFrac * p.frisch);

		int k = model->proddist.size();
		model->matrices->prodmarkov = vector2eigenm(HankUtilities::read_matrix(markov_loc), k, k);
		fix_rounding(model->matrices->prodmarkov);

		for (auto el : model->logprodgrid)
			model->prodgrid.push_back(exp(el));

		model->nprod = model->prodgrid.size();

		// Normalize mean productivity
		double lmean = EigenFunctions::dot(model->prodgrid, model->proddist);

		for (auto& el : model->prodgrid)
			el *= p.meanlabeff / lmean;
	}

	void create_combined_variables(Model* model, const Parameters& p) {
		model->ny = model->nprod * p.nocc;

		model->matrices->ymarkov = MatrixXr::Zero(model->ny, model->ny);
		model->matrices->ymarkovdiag = MatrixXr::Zero(model->ny, model->ny);
		model->yprodgrid.resize(model->ny);
		model->yoccgrid.resize(model->ny);
		model->ydist.resize(model->ny);

		int iy = 0;
		for (int io=0; io<p.nocc; ++io) {
			for (int ip=0; ip<model->nprod; ++ip) {
				model->yprodgrid[iy] = model->prodgrid[ip];
				model->ydist[iy] = model->proddist[ip] * model->occdist[io];
				// yoccgrid[iy] = occgrid[iy];

				for (int ip2=0; ip2<model->nprod; ++ip2) {
					model->matrices->ymarkov(iy, ip2 + model->nprod * io) = model->matrices->prodmarkov(ip, ip2);
				}

				++iy;
			}
		}

		for (int iy=0; iy<model->ny; ++iy)
			model->matrices->ymarkovdiag(iy,iy) = model->matrices->ymarkov(iy,iy);
		
		model->matrices->ymarkovoff = model->matrices->ymarkov - model->matrices->ymarkovdiag;
		model->profsharegrid = model->yprodgrid;
		for (auto& el : model->profsharegrid)
			el /= p.meanlabeff;
	}

	void check_nbl(double bgrid0, const Parameters& p) {
		if ( p.borrowing ) {
			double nbl = -p.lumptransfer / (p.rborr + p.perfectAnnuityMarkets * p.deathrate);

			if (bgrid0 < nbl) {
				throw "Natural borrowing limit violated";
			}
		}
	}

	void fix_rounding(MatrixXr& mat) {
		for (int i=0; i<mat.rows(); ++i) {
			double rowsum = mat.row(i).sum();
			mat(i,i) -= rowsum;
			assert( abs(mat.row(i).sum()) < 1e-7 );
		}
	}

	std::vector<hank_float_type> compute_grid_deltas(
		const std::vector<hank_float_type>& grid, const std::vector<hank_float_type>& dgrid)
	{
		int n = grid.size();
		std::vector<hank_float_type> deltas(n);

		deltas[0] = 0.5 * dgrid[0];
		for (int id=1; id<n-1; ++id)
			deltas[id] = 0.5 * (dgrid[id-1] + dgrid[id]);

		deltas[n-1] = 0.5 * dgrid[n-2];

		return deltas;
	}

	void powerSpacedGrid(double low, double high, double curv, std::vector<hank_float_type>& grid) {
		int n = grid.size();
		HankNumerics::linspace(0.0, 1.0, n, grid);

		for (int i=0; i<n; ++i)
			grid[i] = low + (high - low) * pow(grid[i], 1.0 / curv);
	}

	void adjustPowerSpacedGrid(std::vector<hank_float_type>& grid) {
		if (grid.size() > 10)
			for (int i=0; i<9; ++i)
				grid[i] = i * grid[9] / (10.0 - 1.0);
	}

	void print_value(const std::string& pname, double value, bool insert_endline) {
		std::cout << "  " << pname << " = " << value;

		if ( insert_endline )
			std::cout << '\n';
	}

	void print_value(const std::string& pname, double value) {
		print_value(pname, value, true);
	}

	void check_adjcosts(const Parameters& p, const std::shared_ptr<AdjustmentCosts>& adjcosts) {
		assert( p.kappa_w_fc == adjcosts->kappa_w_fc );
		assert( p.kappa_d_fc == adjcosts->kappa_d_fc );

		for (int i=0; i<5; ++i) {
			assert( p.kappa_w[i] == adjcosts->kappa_w[i] );
			assert( p.kappa_d[i] == adjcosts->kappa_d[i] );
		}
	}
}