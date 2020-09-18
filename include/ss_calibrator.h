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

		int nvals() const {return obj_functions.size();}

};

// void set_from_x(Parameters *p, Model *model, SteadyState *iss) {

// }



#endif