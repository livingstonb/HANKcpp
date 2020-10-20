#ifndef _MODEL_H
#define _MODEL_H

#include <hank_config.h>
#include <string>
#include <vector>
#include <adjustment_costs.h>
#include <parameters.h>
#include <memory>

// Class to store Eigen::Matrix variables
struct ModelMatrices;

// Binds the attributes constructed in Model to const references
class Model {
	public:
		Model(const Parameters& p);

		~Model();

		std::vector<hank_float_type> bgrid, dbgrid, bdelta, agrid, dagrid, adelta, abdelta;

		std::vector<hank_float_type> occgrid, occdist, prodgrid, proddist, logprodgrid, profsharegrid;

		std::vector<hank_float_type> yprodgrid, yoccgrid, ydist, occYsharegrid, occNsharegrid;

		ModelMatrices* matrices = nullptr;

		int nb, na, nocc, nprod, ny, naby, ntot;

		std::vector<int> dims;

		double prodmarkovscale = 1.0;

		std::shared_ptr<AdjustmentCosts> adjcosts = nullptr;

		void make_asset_grids(const Parameters& p);

		void make_occupation_grids(const Parameters& p);

		void create_income_process(const Parameters& p);

		void create_combined_variables(const Parameters& p);

		void check_nbl(const Parameters& p) const;

		std::vector<hank_float_type> get_rb_effective() const;

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

		const Parameters& p;
};

#endif