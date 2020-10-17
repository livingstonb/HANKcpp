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
		ModelBase(const Parameters& p);

		VectorXr bgrid_, dbgrid_, bdelta_;
		VectorXr agrid_, dagrid_, adelta_;
		VectorXr abdelta_;
		VectorXr occgrid_, occdist_;
		VectorXr logprodgrid_, prodgrid_;
		VectorXr proddist_, profsharegrid_;
		VectorXr yprodgrid_, yoccgrid_, ydist_;
		MatrixXr prodmarkov_, ymarkov_, ymarkovdiag_, ymarkovoff_;

		VectorXr occYsharegrid_, occNsharegrid_;

		int nocc_, nprod_;

		AdjustmentCosts adjcosts_;

		void make_asset_grids(const Parameters& p);

		void make_occupation_grids(const Parameters& p);

		void create_income_process(const Parameters& p);

		void create_combined_variables(const Parameters& p);

		void check_nbl(const Parameters& p) const;
};

// Binds the attributes constructed in ModelBase to const references
class Model : private ModelBase {
	public:
		Model(const Parameters& p_);

		const Parameters p;
		const AdjustmentCosts& adjcosts = adjcosts_;
		const VectorXr& bgrid = bgrid_;
		const VectorXr& dbgrid = dbgrid_;
		const VectorXr& bdelta = bdelta_;
		const VectorXr& agrid = agrid_;
		const VectorXr& dagrid = dagrid_;
		const VectorXr& adelta = adelta_;
		const VectorXr& abdelta = abdelta_;
		const VectorXr& occgrid = occgrid_;
		const VectorXr& occdist = occdist_;
		const VectorXr& logprodgrid = logprodgrid_;
		const VectorXr& prodgrid = prodgrid_;
		const VectorXr& proddist = proddist_;
		const VectorXr& profsharegrid = profsharegrid_;
		const MatrixXr& prodmarkov = prodmarkov_;
		const VectorXr& yprodgrid = yprodgrid_;
		const VectorXr& yoccgrid = yoccgrid_;
		const VectorXr& ydist = ydist_;
		const MatrixXr& ymarkov = ymarkov_;
		const MatrixXr& ymarkovdiag = ymarkovdiag_;
		const MatrixXr& ymarkovoff = ymarkovoff_;
		const int nb = p.nb;
		const int na = p.na;
		const int nocc = nocc_;
		const int nprod = nprod_;
		const int ny = nprod_ * nocc_;
		const int naby = p.nb * p.na * nocc_ * nprod_;
		const int ntot = p.nb * p.na * nocc_ * nprod_;
		const std::vector<int> dims = {p.na, p.nb, ny};
		const double prodmarkovscale = 1.0;

		const VectorXr& occYsharegrid = occYsharegrid_;
		const VectorXr& occNsharegrid = occNsharegrid_;

		VectorXr get_rb_effective() const;

		hank_float_type util(hank_float_type c, hank_float_type riskaver) const;

		hank_float_type util1(hank_float_type c, hank_float_type riskaver) const;

		hank_float_type util1inv(hank_float_type u, hank_float_type riskaver) const;

		hank_float_type labdisutil(hank_float_type h, hank_float_type chi) const;

		hank_float_type labdisutil1(hank_float_type h, hank_float_type chi) const;

		hank_float_type labdisutil1inv(hank_float_type du, hank_float_type chi) const;

		double util1BC(double h, double riskaver, double chi, double bdrift, double netwage, double wagescale) const;

		hank_float_type capadjcost1inv(hank_float_type x) const;

		hank_float_type capadjcost1(hank_float_type x) const;

		hank_float_type capadjcost(hank_float_type x) const;

		void print_values() const;

		void assertions() const;
};

#endif