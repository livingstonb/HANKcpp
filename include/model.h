#ifndef _MODEL_H
#define _MODEL_H

#include <hank_config.h>
#include <string>
#include <vector>

#include <memory>

// Class to store Eigen::Matrix variables
struct ModelMatrices;

class Parameters;

class AdjustmentCosts;

class Model {
	private:
		const Parameters& p;

	public:
		Model(const Parameters& p);

		std::vector<hank_float_type> bgrid, dbgrid, bdelta, agrid, dagrid, adelta, abdelta;

		std::vector<hank_float_type> occgrid, occdist, prodgrid, proddist, logprodgrid, profsharegrid;

		std::vector<hank_float_type> yprodgrid, yoccgrid, ydist, occYsharegrid, occNsharegrid;

		std::shared_ptr<ModelMatrices> matrices = nullptr;

		int nb, na, nocc, nprod, ny, naby, ntot;

		std::vector<int> dims;

		double prodmarkovscale = 1.0;

		std::shared_ptr<AdjustmentCosts> adjcosts = nullptr;

		std::vector<hank_float_type> get_rb_effective(hank_float_type rb, hank_float_type rborr) const;

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