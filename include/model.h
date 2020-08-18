#ifndef _MODEL_H
#define _MODEL_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <boost/algorithm/string.hpp>

#include <parameters.h>
#include <utilities.h>
#include <functions.h>
#include <hank.h>
#include <adjustment_costs.h>

class ModelBase
{
	public:
		ModelBase(const Parameters p, const std::string& income_dir) {
			make_asset_grids(p);
			make_occupation_grids(p);
			create_income_process(income_dir, p);
			create_combined_variables(p);

			check_nbl(p);

			adjcosts_ = AdjustmentCosts(p.adjCostRatioMode, p.exponential_adjcosts,
				p.kappa_w_fc, p.kappa_d_fc, p.kappa_w, p.kappa_d);
		}

		double_vector bgrid_;
		double_vector dbgrid_;
		double_vector bdelta_;
		double_vector agrid_;
		double_vector dagrid_;
		double_vector adelta_;
		double_vector occgrid_;
		double_vector occdist_;

		double_vector logprodgrid_;
		double_vector prodgrid_;
		double_vector proddist_;
		double_vector profsharegrid_;
		double_matrix prodmarkov_;

		double_vector yprodgrid_;
		double_vector yoccgrid_;
		double_vector ydist_;
		double_matrix ymarkov_;

		int nocc_;
		int nprod_;

		AdjustmentCosts adjcosts_;

		void make_asset_grids(const Parameters& p);
		void make_occupation_grids(const Parameters& p);
		void create_income_process(const std::string& income_dir, const Parameters& p);
		void create_combined_variables(const Parameters& p);
		void check_nbl(const Parameters& p) const;
};

class Model : private ModelBase {
	public:
		Model(const Parameters p_, const std::string& income_dir)
			: ModelBase(p_, income_dir), p(p_), dims({p_.na, p_.nb, p_.ny}) {};

		const Parameters p;
		const AdjustmentCosts& adjcosts = adjcosts_;

		const double_vector& bgrid = bgrid_;
		const double_vector& dbgrid = dbgrid_;
		const double_vector& bdelta = bdelta_;
		const double_vector& agrid = agrid_;
		const double_vector& dagrid = dagrid_;
		const double_vector& adelta = adelta_;
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

		const int nb = p.nb;
		const int na = p.na;
		const int nocc = nocc_;
		const int nprod = nprod_;
		const int ntot = p.nb * p.na * nocc_ * nprod_;
		const boost_array_shape<double, 3> dims;

		double_vector get_rb_effective() const;
		double util(double c) const;
		double util1(double c) const;
		double util1inv(double u) const;
		double labdisutil(double h, double chi) const;
		double labdisutil1(double u, double chi) const;
		double labdisutil1inv(double du, double chi) const;
		double util1BC(double h, double chi, double bdrift, double netwage, double wagescale) const;
};

#endif