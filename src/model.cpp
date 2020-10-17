#include <model.h>
#include <utilities.h>
#include <hank_numerics.h>
#include <functions.h>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <assert.h>
#include <hank_macros.h>
#include <adjustment_costs.h>

namespace {
	void fix_rounding(MatrixXr& mat);

	VectorXr compute_grid_deltas(const VectorXr& grid, const VectorXr& dgrid);

	void powerSpacedGrid(double low, double high, double curv, VectorXr& grid);

	void adjustPowerSpacedGrid(VectorXr& grid);

	void print_value(const std::string& pname, double value, bool insert_endline);

	void print_value(const std::string& pname, double value);

	void check_adjcosts(const Parameters& p, const AdjustmentCosts& adjcosts);
}

ModelBase::ModelBase(const Parameters& p) {
	make_asset_grids(p);
	make_occupation_grids(p);
	create_income_process(p);
	create_combined_variables(p);

	check_nbl(p);

	double adjcost1max = 1.0e30;

	AdjustmentCosts adjcosts_initial(p.adjCostRatioMode, p.exponential_adjcosts,
		p.kappa_w_fc, p.kappa_d_fc, p.kappa_w, p.kappa_d, adjcost1max, p.dmax);

	adjcost1max = adjcosts_initial.cost1(p.dmax, 1.0);
	adjcosts_ = std::move(AdjustmentCosts(p.adjCostRatioMode, p.exponential_adjcosts,
		p.kappa_w_fc, p.kappa_d_fc, p.kappa_w, p.kappa_d, adjcost1max, p.dmax));
}

void ModelBase::make_asset_grids(const Parameters& p) {
	// Liquid asset
	bgrid_ = VectorXr(p.nb);
	if ( !p.borrowing ) {
		powerSpacedGrid(p.bmin, p.bmax, p.bcurv, bgrid_);
	}
	else {
		VectorXr bgridpos(p.nb_pos);
		powerSpacedGrid(0.0, p.bmax, p.bcurv, bgridpos);

		double nbl = -p.lumptransfer / (p.rborr + p.perfectAnnuityMarkets * p.deathrate);
		double abl = fmax(nbl + p.cmin, p.blim);

		int nn = static_cast<int>(floor(p.nb_neg / 2.0) + 1);
		VectorXr bgridneg(nn);
		powerSpacedGrid(abl, (abl+bgridpos[0])/2.0, p.bcurv_neg, bgridneg);

		bgrid_(seq(0, nn-1)) = bgridneg;
		bgrid_(seq(p.nb_neg, p.nb-1)) = bgridpos;
		for (int i=nn; i<p.nb_neg; ++i)
			bgrid_[i] = bgrid_[p.nb_neg] - (bgrid_[p.nb_neg-i] - bgrid_[0]);
	}
	

	dbgrid_ = bgrid_(seq(1,p.nb-1)) - bgrid_(seq(0,p.nb-2));
	bdelta_ = compute_grid_deltas(bgrid_, dbgrid_);

	// Illiquid asset
	if ( p.oneAssetNoCapital ) {
		agrid_ = VectorXr::Constant(p.na, 0.0);
		dagrid_ = VectorXr::Constant(p.na-1, 1.0);
		adelta_ = VectorXr::Constant(p.na, 1.0);
	}
	else {
		agrid_ = VectorXr(p.na);
		powerSpacedGrid(p.amin, p.amax, p.acurv, agrid_);
		adjustPowerSpacedGrid(agrid_);

		dagrid_ = agrid_(seq(1,p.na-1)) - agrid_(seq(0,p.na-2));
		adelta_ = compute_grid_deltas(agrid_, dagrid_);
	}

	abdelta_.resize(p.nab);
	for (int ia=0; ia<p.na; ++ia)
		for (int ib=0; ib<p.nb; ++ib)
			abdelta_(TO_INDEX_1D(ia, ib, p.na, p.nb)) = adelta_(ia) * bdelta_(ib);
}

void ModelBase::make_occupation_grids(const Parameters& p) {
	// Occupation types
	occYsharegrid_.resize(p.nocc);
	occNsharegrid_.resize(p.nocc);
	occdist_.resize(p.nocc);
	if ( p.nocc == 1 ) {
		occYsharegrid_ << 1;
		occNsharegrid_ << 1;
		occdist_ << 1;
	}
	else if ( p.nocc == 4 ) {
		occYsharegrid_ << 0.325, 0.275, 0.225, 0.175;
		occNsharegrid_ << 0.1, 0.15, 0.25, 0.5;
		occdist_ << 0.25, 0.25, 0.25, 0.25;
	}
	else {
		std::cerr << "Invalid number of occupation points\n";
		throw 0;
	}
	nocc_ = p.nocc;
}

void ModelBase::create_income_process(const Parameters& p) {

	std::string grid_loc = "../input/" + p.income_dir + "/ygrid_combined.txt";
	logprodgrid_ = vector2eigenv(HankUtilities::read_matrix(grid_loc));

	std::string dist_loc = "../input/" + p.income_dir + "/ydist_combined.txt";
	proddist_ = vector2eigenv(HankUtilities::read_matrix(dist_loc));

	std::string markov_loc = "../input/" + p.income_dir + "/ymarkov_combined.txt";

	if ( p.adjustProdGridFrisch )
		logprodgrid_ = logprodgrid_ / (1.0 + p.adjFrischGridFrac * p.frisch);

	int k = proddist_.size();
	prodmarkov_ = vector2eigenm(HankUtilities::read_matrix(markov_loc), k, k);
	fix_rounding(prodmarkov_);

	prodgrid_ = logprodgrid_.array().exp();
	nprod_ = prodgrid_.size();

	// Normalize mean productivity
	double lmean = prodgrid_.dot(proddist_);
	prodgrid_ *= p.meanlabeff / lmean;
}

void ModelBase::create_combined_variables(const Parameters& p) {
	int ny = nprod_ * nocc_;

	ymarkov_ = MatrixXr::Zero(ny, ny);
	ymarkovdiag_ = MatrixXr::Zero(ny, ny);
	yprodgrid_ = VectorXr(ny);
	yoccgrid_ = VectorXr::Zero(ny);
	ydist_ = VectorXr(ny);

	int iy = 0;
	for (int io=0; io<nocc_; ++io) {
		for (int ip=0; ip<nprod_; ++ip) {
			yprodgrid_(iy) = prodgrid_(ip);
			// yoccgrid_(iy) = occgrid_(io);
			ydist_(iy) = proddist_(ip) * occdist_(io);

			for (int ip2=0; ip2<nprod_; ++ip2) {
				ymarkov_(iy, ip2 + nprod_ * io) = prodmarkov_(ip, ip2);
			}

			++iy;
		}
	}

	for (int iy=0; iy<ny; ++iy)
		ymarkovdiag_(iy,iy) = ymarkov_(iy,iy);
	
	ymarkovoff_ = ymarkov_ - ymarkovdiag_;
	profsharegrid_ = yprodgrid_.array() / p.meanlabeff;
}

void ModelBase::check_nbl(const Parameters& p) const {
	if ( p.borrowing ) {
		double nbl = -p.lumptransfer / (p.rborr + p.perfectAnnuityMarkets * p.deathrate);

		if (bgrid_(0) < nbl) {
			throw "Natural borrowing limit violated";
		}
	}
}

Model::Model(const Parameters& p_)
	: ModelBase(p_), p(p_) {
	assertions();

	if ( global_hank_options->print_diagnostics )
		print_values();

	check_adjcosts(p_, adjcosts);
}

VectorXr Model::get_rb_effective() const {
	VectorXr rb_effective, bvec = bgrid;
	rb_effective = bvec.unaryExpr([this](hank_float_type x) -> hank_float_type {
			return (x >= static_cast<hank_float_type>(0.0)) ? p.rb : p.rborr;
		});
	rb_effective = rb_effective.array() + p.perfectAnnuityMarkets * p.deathrate;

	return rb_effective;
}

double Model::util(double c, double riskaver) const {
	return HankFunctions::utility(c, p.prefshock, riskaver);
}

double Model::util1(double c, double riskaver) const {
	return HankFunctions::utility1(c, p.prefshock, riskaver);
}

double Model::util1inv(double u, double riskaver) const {
	return HankFunctions::utility1inv(u, p.prefshock, riskaver);
}

double Model::labdisutil(double h, double chi) const {
	return HankFunctions::labor_disutility(h, p.frisch, chi);
}

double Model::labdisutil1(double h, double chi) const {
	return HankFunctions::labor_disutility1(h, p.frisch, chi);
}

double Model::labdisutil1inv(double du, double chi) const {
	return HankFunctions::labor_disutility1inv(du, p.frisch, chi);
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
		print_value("prodgrid[io]", prodgrid(io));

	for (int io=0; io<nprod; ++io)
		print_value("proddist[io]", proddist(io));

	print_value("E[prod]", prodgrid.dot(proddist));
	
	HankUtilities::horzline();
}

void Model::assertions() const {
	if ( abs(1.0-proddist.sum()) > 1.0e-7 ) {
		std::cerr << "Proddist does not sum to one\n";
		throw 0;
	}

	if ( abs(1.0-occdist.sum()) > 1.0e-7 ) {
		std::cerr << "Occdist does not sum to one\n";
		throw 0;
	}

	if ( abs(1.0-ydist.sum()) > 1.0e-7 ) {
		std::cerr << "ydist does not sum to one\n";
		throw 0;
	}

	ArrayXr rowsums = ymarkov.rowwise().sum();
	if ( rowsums.abs().maxCoeff() > 1.0e-7 ) {
		std::cerr << "Markov ytrans matrix rows do not sum to zero\n";
		throw 0;
	}
}

namespace {
	void fix_rounding(MatrixXr& mat) {
		for (int i=0; i<mat.rows(); ++i) {
			double rowsum = mat.row(i).sum();
			mat(i,i) -= rowsum;
			assert( abs(mat.row(i).sum()) < 1e-7 );
		}
	}

	VectorXr compute_grid_deltas(const VectorXr& grid, const VectorXr& dgrid) {
		int n = grid.size();
		int n_d = n - 1;
		VectorXr deltas(n);

		deltas(0) = 0.5 * dgrid(0);
		deltas(seq(1,n-2)) = 0.5 * (dgrid(seq(0,n_d-2)) + dgrid(seq(1,n_d-1)));
		deltas(n-1) = 0.5 * dgrid(n_d-1);

		return deltas;
	}

	void powerSpacedGrid(double low, double high, double curv, VectorXr& grid) {
		int n = grid.size();
		HankNumerics::linspace(0.0, 1.0, n, grid);

		for (int i=0; i<n; ++i)
			grid[i] = low + (high - low) * pow(grid[i], 1.0 / curv);
	}

	void adjustPowerSpacedGrid(VectorXr& grid) {
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

	void check_adjcosts(const Parameters& p, const AdjustmentCosts& adjcosts) {
		assert( p.kappa_w_fc == adjcosts.kappa_w_fc );
		assert( p.kappa_d_fc == adjcosts.kappa_d_fc );

		for (int i=0; i<5; ++i) {
			assert( p.kappa_w[i] == adjcosts.kappa_w[i] );
			assert( p.kappa_d[i] == adjcosts.kappa_d[i] );
		}
	}
}