#include <model.h>

namespace {
	void fix_rounding(double_matrix& mat)
	{
		for (int i=0; i<mat.rows(); ++i)
			mat(i,i) = mat(i,i) - mat.row(i).sum();
	}

	double_vector compute_grid_deltas(const double_vector& grid, const double_vector& dgrid)
	{
		int n = grid.size();
		double_vector deltas(n);

		deltas(0) = 0.5 * dgrid(0);
		deltas(seq(1,n-2)) = 0.5 * (dgrid(seq(0,last-1)) + dgrid(seq(1,last)));
		deltas(n-1) = 0.5 * dgrid(last);

		return deltas;
	}
}

void ModelBase::make_asset_grids(const Parameters& p) {
	// Liquid asset
	bgrid_ = double_vector(p.nb);
	if ( !p.borrowing ) {
		powerSpacedGrid(p.bmin, p.bmax, p.bcurv, bgrid_);
	}
	else {
		double_vector bgridpos(p.nb_pos);
		powerSpacedGrid(0.0, p.bmax, p.bcurv, bgridpos);

		double nbl = -p.lumptransfer / (p.rborr + p.perfectAnnuityMarkets * p.deathrate);
		double abl = fmax(nbl + p.cmin, p.blim);

		int nn = static_cast<int>(floor(p.nb_neg / 2.0) + 1);
		double_vector bgridneg(nn);
		powerSpacedGrid(abl, (abl+bgridpos[0])/2.0, p.bcurv, bgridneg);

		bgrid_(seq(0, nn-1)) = bgridneg;
		bgrid_(seq(p.nb_neg, p.nb-1)) = bgridpos;
		for (int i=nn; i<p.nb_neg; ++i)
			bgrid_[i] = bgrid_[p.nb_neg] - (bgrid_[p.nb_neg+2-i] - bgrid_[0]);
	}
	

	dbgrid_ = bgrid_(seq(1,p.nb-1)) - bgrid_(seq(0,p.nb-2));
	bdelta_ = compute_grid_deltas(bgrid_, dbgrid_);

	// Illiquid asset
	agrid_ = double_vector(p.na);
	powerSpacedGrid(p.amin, p.amax, p.acurv, agrid_);
	adjustPowerSpacedGrid(agrid_);

	dagrid_ = agrid_(seq(1,p.na-1)) - agrid_(seq(0,p.na-2));
	adelta_ = compute_grid_deltas(agrid_, dagrid_);
}

void ModelBase::make_occupation_grids(const Parameters& p) {
	occgrid_ = double_vector(p.nocc);
	occdist_ = double_vector(p.nocc);
	nocc_ = p.nocc;

	if (p.nocc == 1) {
		std::fill(occgrid_.begin(), occgrid_.end(), 0.0);
		std::fill(occdist_.begin(), occdist_.end(), 1.0);
	}
	else {
		// S_N / (S_N + S_Y)
		double lshareNY = (1.0 - p.alpha_N) * p.drs_N
			/ ((p.elast - 1.0) * (1.0 - p.alpha_Y) * p.drs_Y + (1.0 - p.alpha_N) * p.drs_N);

		if (lshareNY == 0.0) {
			// No labor income accrues to N-type
			std::fill(occgrid_.begin(), occgrid_.end(), 0.0);
			std::fill(occdist_.begin(), occdist_.end(), 1.0 / p.nocc);
		}
		else if (lshareNY == 1.0) {
			// All labor income accrues to N-type
			std::fill(occgrid_.begin(), occgrid_.end(), 1.0);
			std::fill(occdist_.begin(), occdist_.end(), 1.0 / p.nocc);
		}
		else {
			// Equally spaced in [0, 1], midpoints of intervals
			double lwidth = 1.0 / p.nocc;
			occgrid_[0] = 0.5 * lwidth;
			occgrid_[p.nocc-1] = 1.0 - 0.5 * lwidth;
			if (p.nocc >= 2) {
				for (int i=1; i<p.nocc-1; ++i) {
					occgrid_[i] = occgrid_[i-1] + lwidth;
				}
			}

			// Distribution has CDF x ^ par. Choose par to target mean a
			double lWNtoWY = 1.5;
			double lmeanocc = lshareNY / (lshareNY + (1.0 - lshareNY) * lWNtoWY);
			double lpar = lmeanocc / (1.0 - lmeanocc);
			for (int i=0; i<p.nocc; ++i) {
				occdist_[i] = pow(occgrid_[i] + 0.5 * lwidth, lpar)
					- pow(occgrid_[i] - 0.5 * lwidth, lpar);
			}
		}
	}
}

void ModelBase::create_income_process(
	const std::string& income_dir, const Parameters& p) {

	std::string grid_loc = "input/" + income_dir + "/ygrid_combined.txt";
	logprodgrid_ = vector2eigenv(read_matrix(grid_loc));

	std::string dist_loc = "input/" + income_dir + "/ydist_combined.txt";
	proddist_ = vector2eigenv(read_matrix(dist_loc));

	std::string markov_loc = "input/" + income_dir + "/ymarkov_combined.txt";
	int k = proddist_.size();
	prodmarkov_ = vector2eigenm(read_matrix(markov_loc), k, k);
	fix_rounding(prodmarkov_);

	prodgrid_ = logprodgrid_.array().exp();
	nprod_ = prodgrid_.size();

	// Normalize mean productivity
	double lmean = prodgrid_.dot(proddist_);
	prodgrid_ = p.meanlabeff * prodgrid_ / lmean;
}

void ModelBase::create_combined_variables(const Parameters& p) {
	int iy, io, ip, iy2, io2, ip2;
	int ny = nprod_ * nocc_;

	double_vector occfromy(ny);
	double_vector prodfromy(ny);
	double_matrix yfromoccprod(nocc_, nprod_);

	iy = 0;
	for (int io=0; io<nocc_; ++io) {
		for (int ip=0; ip<nprod_; ++ip) {
			occfromy(iy) = io;
			prodfromy(iy) = ip;
			yfromoccprod(io,ip) = iy;
			++iy;
		}
	}

	ymarkov_ = double_matrix::Zero(ny, ny);
	ymarkovdiag_ = double_matrix::Zero(ny, ny);
	yprodgrid_ = double_vector(ny);
	yoccgrid_ = double_vector(ny);
	ydist_ = double_vector(ny);

	for (int iy=0; iy<ny; ++iy) {
		io = occfromy(iy);
		ip = prodfromy(iy);
		yprodgrid_(iy) = prodgrid_(ip);
		yoccgrid_(iy) = occgrid_(io);
		ydist_(iy) = proddist_(ip) * occdist_(io);

		for (iy2=0; iy2<ny; ++iy2) {
			io2 = occfromy(iy2);
			ip2 = prodfromy(iy2);

			if (io == io2)
				ymarkov_(iy,iy2) = prodmarkov_(ip,ip2);
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

double_vector Model::get_rb_effective() const
{
	double_vector rb_effective = bgrid;
	rb_effective = rb_effective.unaryExpr([this](double x) {
			return (x >= 0.0) ? p.rb : p.rborr;
		});
	rb_effective = rb_effective.array() + p.perfectAnnuityMarkets * p.deathrate;

	return rb_effective;
}

double Model::util(double c) const {
	return HankFunctions::utility(c, p.prefshock, p.riskaver);
}

double Model::util1(double c) const {
	return HankFunctions::utility1(c, p.prefshock, p.riskaver);
}

double Model::util1inv(double u) const {
	return HankFunctions::utility1inv(u, p.prefshock, p.riskaver);
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

double Model::util1BC(double h, double chi, double bdrift, double netwage, double wagescale) const {
	return labdisutil1(h, chi) - util1(bdrift + h * netwage) * netwage * wagescale * p.labwedge;
}