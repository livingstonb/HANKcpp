#include <model.h>

void ModelBase::make_grids(const Parameters& p) {
	// Liquid asset
	bgrid_ = double_vector(p.nb);
	powerSpacedGrid(p.nb, p.bmin, p.bmax, p.bcurv, bgrid_);

	bdelta_ = compute_grid_deltas(bgrid_);

	// Illiquid asset
	agrid_ = double_vector(p.na);
	powerSpacedGrid(p.na, p.amin, p.amax, p.acurv, agrid_);
	adjustPowerSpacedGrid(agrid_);

	adelta_ = compute_grid_deltas(agrid_);

	// Occupations
	auto occgrids = occupationGrid(p);
	occgrid_ = vector2eigenv(occgrids.first);
	occdist_ = vector2eigenv(occgrids.second);
	nocc_ = occgrid_.size();
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

void ModelBase::create_combined_variables() {
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
}

double_vector Model::get_rb_effective() const
{
	double_vector rb_effective = bgrid;
	rb_effective.unaryExpr([=, *this](double x) {
			return (x >= 0.0) ? rb : rborr;
		});
	rb_effective = rb_effective.array() + perfectAnnuityMarkets * deathrate;

	return rb_effective;
}

std::vector<double> read_matrix(const std::string& file_loc)
{
	std::string line, word;
	std::ifstream yfile;
	std::size_t current, previous;
	yfile.open(file_loc.data(), std::ios::in);

	std::vector<double> out;

	while ( getline(yfile, line) ) {
		previous = 0;
    	current = find_multiple(line, 0);
    	if (current != std::string::npos) {
	    	while (current != std::string::npos) {
	    		if (current > 0) {
		    		word = line.substr(previous, current - previous);
		    		out.push_back(std::stod(word));
			    }

			    previous = current + 1;
		    	current = find_multiple(line, previous);
		    }

		    word = line.substr(previous, current - previous);
		    out.push_back(std::stod(word));
		}
		else {
			out.push_back(std::stod(line));
		}
	}
	yfile.close();

	return out;
}

std::size_t find_multiple(const std::string& line, int pos)
{
	std::size_t t1, t2;

	t1 = line.find("  ", pos);
	t2 = line.find(" -", pos);

	if ((t1 < t2) & (t1 != std::string::npos)) {
		return t1;
	}
	else {
		return t2;
	}
}

void fix_rounding(double_matrix& mat)
{
	for (int i=0; i<mat.rows(); ++i)
		mat(i,i) = mat(i,i) - mat.row(i).sum();
}

double_vector compute_grid_deltas(const double_vector& grid)
{
	int n = grid.size();
	double_vector dgrid = grid(seq(1,last)) - grid(seq(0,last-1));
	double_vector deltas(n);

	deltas(0) = 0.5 * dgrid(0);
	deltas(seq(1,n-2)) = 0.5 * (dgrid(seq(0,last-1)) + dgrid(seq(1,last)));
	deltas(n-1) = 0.5 * dgrid(last);

	return deltas;
}