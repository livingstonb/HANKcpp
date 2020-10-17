#ifndef _SS_DEVIATIONS_H
#define _SS_DEVIATIONS_H

#include <hank_config.h>
#include <hank_types.h>
#include <parameters.h>
#include <model.h>
#include <equilibrium.h>
#include <distribution_statistics.h>
#include <vector>
#include <string>

namespace HANKCalibration {

struct SSCalibrationArgs {
	SSCalibrationArgs() {}
	
	SSCalibrationArgs(const Parameters *p_, const Model *model_,
		const DistributionStatistics *stats_, const EquilibriumElement *iss_) {
		p = p_;
		model = model_;
		stats = stats_;
		iss = iss_;
	}

	const Parameters *p = NULL;

	const Model *model = NULL;

	const DistributionStatistics *stats = NULL;

	const EquilibriumElement *iss = NULL;
};

using deviation_fn_type = std::function<double(const SSCalibrationArgs&)>;

class SSCalibrator {
	public:
		SSCalibrator() {}

		void setup(const Parameters &p);

		std::vector<deviation_fn_type> obj_functions;

		void fill_fvec(const SSCalibrationArgs& args, hank_float_type fvec[]) const;

		void fill_xguess(const Parameters &p, const Model& model, hank_float_type xvec[]);

		void update_params(Parameters *p, const hank_float_type *xvec) const;

		void update_ss(const Parameters* p, EquilibriumElement *iss, const hank_float_type *xvec) const;\

		void print_fvec(hank_float_type fvec[]) const;

		int nmoments;

		void check_size(int ix) const;

		void perform_calibrator_assertions() const;

		bool calibrateLaborDisutility = true;

		bool calibrateRb = true;

		bool calibrateDiscountRate = true;

		std::vector<std::string> moment_descriptions;

		std::vector<int> ix_labor_occ;

		int ix_rho = -1;

		int ix_rb = -1;

		int ix_chi = -1;

		int ix_capital = -1;
};

using ObjectPointers = UniquePtrContainer<Parameters, Model, EquilibriumElement, DistributionStatistics, SSCalibrator>;

int initial_state_state_obj_fn(void* args_void_ptr, int n, const hank_float_type *x, hank_float_type *fvec, int /* iflag */ );

}

#endif