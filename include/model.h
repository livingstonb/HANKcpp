#ifndef _MODEL_H
#define _MODEL_H

#include <hank_config.h>
#include <string>
#include <vector>
#include <hank_eigen_dense.h>
#include <adjustment_costs.h>
#include <parameters.h>

// Constructs the grids and provides a container with public attribute access
class ModelBase {
	public:
		ModelBase(const Parameters& p, const std::string& income_dir);

		double_vector bgrid_, dbgrid_, bdelta_;
		double_vector agrid_, dagrid_, adelta_;
		double_vector abdelta_;
		double_vector occgrid_, occdist_;
		double_vector logprodgrid_, prodgrid_;
		double_vector proddist_, profsharegrid_;
		double_vector yprodgrid_, yoccgrid_, ydist_;
		double_matrix prodmarkov_, ymarkov_, ymarkovdiag_, ymarkovoff_;

		double_vector occYsharegrid_, occNsharegrid_;

		int nocc_, nprod_;

		AdjustmentCosts adjcosts_;

		void make_asset_grids(const Parameters& p);

		void make_occupation_grids(const Parameters& p);

		void create_income_process(const std::string& income_dir, const Parameters& p);

		void create_combined_variables(const Parameters& p);

		void check_nbl(const Parameters& p) const;
};

// Binds the attributes constructed in ModelBase to const references
class Model : private ModelBase {
	public:
		Model(const Parameters& p_, const std::string& income_dir);

		const Parameters p;
		const AdjustmentCosts& adjcosts = adjcosts_;
		const double_vector& bgrid = bgrid_;
		const double_vector& dbgrid = dbgrid_;
		const double_vector& bdelta = bdelta_;
		const double_vector& agrid = agrid_;
		const double_vector& dagrid = dagrid_;
		const double_vector& adelta = adelta_;
		const double_vector& abdelta = abdelta_;
		const double_vector& occgrid = occgrid_;
		const double_vector& occdist = occdist_;
		const double_vector& logprodgrid = logprodgrid_;
		const double_vector& prodgrid = prodgrid_;
		const double_vector& proddist = proddist_;
		const double_vector& profsharegrid = profsharegrid_;
		const double_matrix& prodmarkov = prodmarkov_;
		const double_vector& yprodgrid = yprodgrid_;
		const double_vector& yoccgrid = yoccgrid_;
		const double_vector& ydist = ydist_;
		const double_matrix& ymarkov = ymarkov_;
		const double_matrix& ymarkovdiag = ymarkovdiag_;
		const double_matrix& ymarkovoff = ymarkovoff_;
		const int nb = p.nb;
		const int na = p.na;
		const int nocc = nocc_;
		const int nprod = nprod_;
		const int ny = nprod_ * nocc_;
		const int naby = p.nb * p.na * nocc_ * nprod_;
		const int ntot = p.nb * p.na * nocc_ * nprod_;
		const std::vector<int> dims = {p.na, p.nb, ny};
		const double prodmarkovscale = 1.0;

		const double_vector& occYsharegrid = occYsharegrid_;
		const double_vector& occNsharegrid = occNsharegrid_;

		double_vector get_rb_effective() const;

		double util(double c) const;

		double util1(double c) const;

		double util1inv(double u) const;

		double labdisutil(double h, double chi) const;

		double labdisutil1(double h, double chi) const;

		double labdisutil1inv(double du, double chi) const;

		double util1BC(double h, double chi, double bdrift, double netwage, double wagescale) const;

		void print_values() const;

		void assertions() const;
};

#endif