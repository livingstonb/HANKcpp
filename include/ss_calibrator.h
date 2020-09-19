#ifndef _SS_DEVIATIONS_H
#define _SS_DEVIATIONS_H

#include <hank_types.h>
#include <parameters.h>
#include <model.h>
#include <steady_state.h>
#include <distribution_statistics.h>
#include <vector>

struct SSCalibrationArgs {
	SSCalibrationArgs() {}
	
	SSCalibrationArgs(const Parameters *p_, const Model *model_,
		const DistributionStatistics *stats_, const SteadyState *iss_) {
		p = p_;
		model = model_;
		stats = stats_;
		iss = iss_;
	}

	const Parameters *p = NULL;
	const Model *model = NULL;
	const DistributionStatistics *stats = NULL;
	const SteadyState *iss = NULL;
};

using deviation_fn_type = std::function<double(const SSCalibrationArgs&)>;

class SSCalibrator {
	public:
		SSCalibrator(const Parameters &p);

		std::vector<deviation_fn_type> obj_functions;

		void fill_fvec(const SSCalibrationArgs& args, double fvec[]) const;

		void fill_xguess(const Parameters &p, const Model& model, double xvec[]);

		void update_params(Parameters *p, double xvec[]) const;

		void update_ss(const Parameters& p, SteadyState *iss, double xvec[]) const;

		int nmoments() const {return obj_functions.size();}

		void check_size(int ix) const;

		std::vector<int> ix_occdist;

		int ix_rho = -1;

		int ix_rb = -1;

		int ix_chi = -1;

		int ix_capital = -1;
};

#endif