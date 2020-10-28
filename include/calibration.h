#ifndef _SS_DEVIATIONS_H
#define _SS_DEVIATIONS_H

#include <hank_config.h>
#include <hank.h>
#include <parameters.h>
#include <model.h>
#include <equilibrium.h>
#include <distribution_statistics.h>
#include <vector>
#include <string>
#include <functional>

namespace HANKCalibration {

class CalibrationArgs;

class SSCalibrator;

using deviation_fn_type = std::function<double(const CalibrationArgs&)>;

using ObjectPointers = HANK::UniquePtrContainer<Parameters, Model, EquilibriumInitial, DistributionStatistics, SSCalibrator>;

class SSCalibrator
{
	public:
		SSCalibrator() {}

		void setup(const Parameters &p);

		std::vector<deviation_fn_type> obj_functions, variable_values;

		void fill_fvec(const CalibrationArgs& args, hank_float_type fvec[]) const;

		void fill_xguess(const Parameters &p, const Model& model, hank_float_type xvec[]);

		void update_params(Parameters *p, const hank_float_type *xvec) const;

		void update_ss(const Parameters* p, EquilibriumInitial *iss, const hank_float_type *xvec) const;\

		void print_fvec(hank_float_type fvec[]) const;

		int nmoments, iter = 0;

		bool printDetailed = true;

		bool calibrateLaborDisutility = true;

		bool calibrateRb = true;

		bool calibrateDiscountRate = true;

		std::vector<std::string> moment_descriptions, variable_names;

		std::vector<int> ix_labor_occ;

		int ix_rho = -1;

		int ix_rb = -1;

		int ix_chi = -1;

		int ix_capital = -1;
};

class CalibrationArgs
{
	public:
		CalibrationArgs(const ObjectPointers& ptrs)
			: p(*ptrs.ptr1), model(*ptrs.ptr2), iss(*ptrs.ptr3), stats(*ptrs.ptr4) {}

		const Parameters& p;

		const Model& model;

		const EquilibriumInitial& iss;

		const DistributionStatistics& stats;
};

int initial_steady_state_obj_fn(void* args_void_ptr, int n, const hank_float_type *x, hank_float_type *fvec, int /* iflag */ );

}

#endif